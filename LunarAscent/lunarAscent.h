/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/Astrodynamics/Aerodynamics/hypersonicLocalInclinationAnalysis.h>
#include <Tudat/Mathematics/GeometricShapes/capsule.h>
#include <Tudat/Mathematics/Statistics/randomVariableGenerator.h>
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
using namespace tudat;

namespace tudat_applications
{
namespace PropagationOptimization2020
{

std::shared_ptr< ThrustAccelerationSettings > getThrustAccelerationModelFromParameters(
        std::vector< double >& thrustParameters,
        const simulation_setup::NamedBodyMap bodyMap,
        const double initialTime,
        const double constantSpecificImpulse );

/*!
 *  Class to compute the thrust direction and magnitude for the lunar ascent vehicle. The current inputs set a
 *  constant thrust magnitude, and a thrust direction in y-(-z) plane of the vertical frame define linearly in time using
 *  equispaced nodes. These settings are to be modified for the assignment.
 */
class LunarAscentThrustGuidance
{
public:

    //! Contructor
    /*!
     * Contructor
     * \param vehicleBody Body object for the ascent vehicle
     * \param initialTime Start time of the propagatiin
     * \param parameterVector Vector of independent variables to be used for thrust parameterization:
     *   - Entry 0: Constant thrust magnitude
     *   - Entry 1: Constant spacing in time between nodes
     *   - Entry 2-6: Thrust angle theta, at five nodes
     */
    LunarAscentThrustGuidance(
            const std::shared_ptr< Body > vehicleBody,
            const double initialTime,
            const std::vector< double > parameterVector ):
        vehicleBody_( vehicleBody ),
        parameterVector_( parameterVector )
    {
        // Retrieve parameters of thrust profile
        thrustMagnitude_ = parameterVector_.at( 0 );
        timeInterval_ = parameterVector_.at( 1 );

        // Create interpolator for thrust angle
        double currentTime = initialTime;
        for( unsigned int i = 0; i < parameterVector_.size( ) - 2; i++ )
        {
            thrustAngleMap_[ currentTime ] = parameterVector_.at( i + 2 );
            currentTime += timeInterval_;
        }
        thrustAngleInterpolator_ = createOneDimensionalInterpolator(
                    thrustAngleMap_, std::make_shared< InterpolatorSettings >(
                        linear_interpolator, huntingAlgorithm, false, use_boundary_value ) );
    }

    //! Function that computes the inertial thrust direction for each state derivative function evaluation
    Eigen::Vector3d getCurrentThrustDirection( const double currentTime )
    {
        // Retrieve thrust angle
        double currentThrustAngle = thrustAngleInterpolator_->interpolate( currentTime );

        // Set thrust in V-frame
        Eigen::Vector3d thrustDirectionInVerticalFrame =
                ( Eigen::Vector3d( ) << 0.0, std::sin( currentThrustAngle ), -std::cos( currentThrustAngle ) ).finished( );

        // Retrieve rotation from V-frame to inertial frame
        Eigen::Quaterniond verticalToInertialFrame =
                vehicleBody_->getFlightConditions( )->getAerodynamicAngleCalculator( )->getRotationQuaternionBetweenFrames(
                    vertical_frame, inertial_frame );

        // Return thrust direction
        return verticalToInertialFrame * thrustDirectionInVerticalFrame;

    }

    //! Function that computes the thrust magnitude for each state derivative function evaluation
    double getCurrentThrustMagnitude( const double currentTime )
    {
        return thrustMagnitude_;
    }

private:

    //! Object containing properties of the vehicle
    std::shared_ptr< Body > vehicleBody_;

    //! Parameter containing the solution parameter vector
    std::vector< double > parameterVector_;

    //! Map containing the thrust (value) as a function of time (key)
    std::map< double, double > thrustAngleMap_;

    //! Object that interpolates the thrust as a function of time
    std::shared_ptr< OneDimensionalInterpolator< double, double > > thrustAngleInterpolator_;

    //! Constant time between thrust angle nodes
    double timeInterval_;

    //! Constant magnitude of the thrust
    double thrustMagnitude_;

};

class LunarAscentProblem
{
public:

    // Constructor that we will actually use
    LunarAscentProblem(
            const simulation_setup::NamedBodyMap bodyMap,
            const std::shared_ptr< IntegratorSettings< > > integratorSettings,
            const std::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings,
            const double initialTime,
            const double constantSpecificImpulse = 300.0 );

    // Standard constructor
    LunarAscentProblem( )
    {
    }

    //! Function to retrieve the map with the propagated state history of the last run
    std::map< double, Eigen::VectorXd > getLastRunPropagatedStateHistory( ) const
    {
        return dynamicsSimulator_->getEquationsOfMotionNumericalSolution( );
    }

    //! Function to retrieve the map with the dependent variable history of the last run
    std::map< double, Eigen::VectorXd > getLastRunDependentVariableHistory( ) const
    {
        return dynamicsSimulator_->getDependentVariableHistory( );
    }

    //! Function to retrieve a shared pointer to the dynamics simulator of the last run
    std::shared_ptr< SingleArcDynamicsSimulator< > > getLastRunDynamicsSimulator( )
    {
        return dynamicsSimulator_;
    }


    //! Fitness function, called to run the simulation. In this form, it is compatible with
    //! the Pagmo optimization library.
    //! \param shapeParameters Decision vector, containing the shape parameters for which
    //! the propagation needs to be run, to find the corresponding fitness value.
    //! \return Returns the vector with doubles, describing the fitness belonging
    //! to the shape parameters passed to the function. Currently returns an empty
    //! vector.
    //!
    std::vector< double > fitness( std::vector< double >& x ) const;



private:

    //! Instance variable holding the body map for the simulation
    mutable simulation_setup::NamedBodyMap bodyMap_;

    //! Object holding the integrator settings
    std::shared_ptr< IntegratorSettings< > > integratorSettings_;

    //! Object holding the propagator settings
    std::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings_;

    //! Object holding the translational state propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalStatePropagatorSettings_;

    //! Instance variable containing the start epoch of the simulation
    double initialTime_;

    //! Variable containing the constant specific impulse in seconds
    double constantSpecificImpulse_;

    //! Object holding the dynamics simulator
    mutable std::shared_ptr<SingleArcDynamicsSimulator< > > dynamicsSimulator_;

    //! Map holding the propagated state history
    mutable std::map< double, Eigen::VectorXd > propagatedStateHistory;

    //! Map holding the dependent variable history
    mutable std::map< double, Eigen::VectorXd > dependentVariableHistory;

};

} // Namespace tudat_applications
} // Namespace PropagationOptimization2020
