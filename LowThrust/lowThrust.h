/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/hodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/baseFunctionsHodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/createBaseFunctionHodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/getRecommendedBaseFunctionsHodographicShaping.h"
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "../applicationOutput.h"

using namespace tudat::ephemerides;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::propagators;
using namespace tudat::aerodynamics;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace tudat::mathematical_constants;
using namespace tudat::reference_frames;
using namespace tudat::shape_based_methods;
using namespace tudat::low_thrust_trajectories;
using namespace tudat;

namespace tudat_applications
{
namespace PropagationOptimization2020
{

//! Function to retrieve time-of-flight from trajectory parameters
double getTrajectoryTimeOfFlight( const std::vector< double >& trajectoryParameters);

//! Function to inital time from trajectory parameters (taking into account a buffer). This buffer is *added* to the
//! departure time from trajectoryParameters
double getTrajectoryInitialTime( const std::vector< double >& trajectoryParameters, const double bufferTime = 0.0 );

//! Function to final time from trajectory parameters (taking into account a buffer). This buffer is *subtracted* to the
//! arrival time from trajectoryParameters
double getTrajectoryFinalTime( const std::vector< double >& trajectoryParameters, const double bufferTime = 0.0 );


//! Get the propagation termination settings for the lunar ascent.
/*!
 * \param simulationStartEpoch Start time of the simulation in seconds.
 * \param maximumDuration Time in seconds, specifying the maximum time duration before which the
 * simulation should stop.
 * \param terminationAltitude Altitude in meters, specifying the maximum altitude before which the
 * simulation should stop.
 * \param vehicleDryMass Dry mass of the vehicle in kg. This is value is used to create a termination
 * condition that mandates the simulation to stop once all fuel has been used up.
 * \return Shared pointer to the PropagationTerminationSettings object.
 */
std::shared_ptr< PropagationTerminationSettings > getPropagationTerminationSettings(
        const std::vector< double >& trajectoryParameters,
        const double targetDistance,
        const double trajectoryTimeBuffer );

//! Function that creates the object that computes the semi-analytical hodographic shaping method
/*!
 *  Function that creates the object that computes the semi-analytical hodographic shaping method
 *  \param trajectoryParameters Parameters that define the trajectory shape (see comments with main function)
 *  \param bodyMap List of body objects
 *  \return Object that computes the hodographic-shape-based trajectory
 */
std::shared_ptr< HodographicShaping > createHodographicShapingObject(
        std::vector< double >& trajectoryParameters,
        const simulation_setup::NamedBodyMap bodyMap );

//! Function that generates thrust acceleration model from thrust parameters
std::shared_ptr< ThrustAccelerationSettings > getThrustAccelerationSettingsFromParameters(
        std::vector< double >& trajectoryParameters,
        const simulation_setup::NamedBodyMap bodyMap );

//! Function to retrieve the semi-analytical calculation of the state along the shape-based trajectory at a given time
Eigen::Vector6d getHodographicLowThrustStateAtEpoch(
        std::vector< double >& trajectoryParameters,
        const simulation_setup::NamedBodyMap bodyMap,
        const double evaluationTime );

class LowThrustProblem
{
public:

    //! Constructor for the problem class
    /*!
     * Constructor for the problem class
     * \param bodyMap List of body objects
     * \param integratorSettings Settings for numerical integrator
     * \param propagatorSettings Settings for propagation of translational state and mass
     * \param specificImpulse Constant specific impulse of thrust
     * \param minimumMarsDistance Minimum permissible Mars distance (for propagation termination)
     * \param timeBuffer Amount of time after shape-based departure time at which to start propagation
     * \param performPropagation Boolean denoting whether to propagate dynamics numerically
     */
    LowThrustProblem(
            const simulation_setup::NamedBodyMap bodyMap,
            const std::shared_ptr< IntegratorSettings< > > integratorSettings,
            const std::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings,
            double specificImpulse,
            double minimumMarsDistance,
            double timeBuffer,
            const bool performPropagation = true );

    //! Default constructor
    LowThrustProblem( ){ }

    //! Function to retrieve the map with the propagated state history computed at last call of fitness function
    std::map< double, Eigen::VectorXd > getLastRunPropagatedStateHistory( ) const
    {
        return dynamicsSimulator_->getEquationsOfMotionNumericalSolution( );
    }

    //! Function to retrieve the map with the dependent variable history computed at last call of fitness function
    std::map< double, Eigen::VectorXd > getLastRunDependentVariableHistory( ) const
    {
        return dynamicsSimulator_->getDependentVariableHistory( );
    }

    //! Function to the dynamics simulator, as created during last call of fitness function
    std::shared_ptr< SingleArcDynamicsSimulator< > > getLastRunDynamicsSimulator( )
    {
        return dynamicsSimulator_;
    }

    //! Function to compute propagate the dynamics of the vehicle defined by the trajectoryParameters
    /*!
     *  Function to compute propagate the dynamics of the vehicle defined by the trajectoryParameters. This function updates
     *  all relevant settings and properties to the new values of these parameters.
     *
     *  NOTE: Presently no fitness is computed, this must be modified during the group assignment
     *
     *  \param thrustParameters Values of parameters defining the trajectory (see main function)
     *  \return Fitness (undefined)
     */
    std::vector< double > fitness( std::vector< double >& trajectoryParameters ) const;

    std::shared_ptr< HodographicShaping > getHodographicShaping( )
    {
        return hodographicShaping_;
    }

private:

    //! Variable holding the body map for the simulation
    mutable simulation_setup::NamedBodyMap bodyMap_;

    //! Object holding the integrator settings
    std::shared_ptr< IntegratorSettings< > > integratorSettings_;

    //! Object holding the propagator settings
    std::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings_;

    //! Object holding the translational state propagator settings (part of the propagatorSettings_)
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalStatePropagatorSettings_;

    //! Constant specific impulse of thrust
    double specificImpulse_;

    //! Minimum permissible Mars distance (for propagation termination)
    double minimumMarsDistance_;

    //! Amount of time after shape-based departure time at which to start propagation
    double timeBuffer_;

    //! Boolean denoting whether to propagate dynamics numerically
    bool performPropagation_;

    //! Object that computes the shape-based trajectory semi-analytically
    mutable std::shared_ptr< HodographicShaping > hodographicShaping_;

    //! Object holding the dynamics simulator, as created during last call of fitness function
    mutable std::shared_ptr<SingleArcDynamicsSimulator< > > dynamicsSimulator_;


};

} // Namespace tudat_applications
} // Namespace PropagationOptimization2020
