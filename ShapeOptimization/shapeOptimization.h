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
using namespace tudat;

namespace tudat_applications
{
namespace PropagationOptimization2020
{

//! Function that creates an aerodynamic database for a capsule, based on a set of shape parameters (see .cpp file for details)
std::shared_ptr< HypersonicLocalInclinationAnalysis > getCapsuleCoefficientInterface(
        const std::shared_ptr< geometric_shapes::Capsule > capsule,
        const std::string directory,
        const std::string filePrefix );

//! Function to set Capsule properties related to vehicle shape (mass, aerodynamic coefficients)
void setVehicleShapeParameters(
        const std::vector< double >& decisionVariables,
        const NamedBodyMap& bodyMap,
        const double vehicleDensity );

//! Function to add the Capsule, and associated shape properties, to the body map.
void addCapsuleToBodyMap( NamedBodyMap& bodyMap,
                          const std::vector< double >& decisionVariables,
                          const double vehicleDensity );

//! Class to set the aerodynamic angles of the capsule (default: all angles 0, angle-of-attack constant at given value)
class CapsuleAerodynamicGuidance: public aerodynamics::AerodynamicGuidance
{
public:

    //! Constructor
    CapsuleAerodynamicGuidance(
            const NamedBodyMap bodyMap_,
            const double fixedAngleOfAttack ):bodyMap_( bodyMap_ ), fixedAngleOfAttack_( fixedAngleOfAttack )
    {

    }

    //! The aerodynamic angles are to be computed here
    void updateGuidance( const double time )
    {
        currentAngleOfAttack_ = fixedAngleOfAttack_;
        currentAngleOfSideslip_ = 0.0;
        currentBankAngle_ = 0.0;

    }

private:

    //! List of body objects that constitute the environment
    NamedBodyMap bodyMap_;

    //! Fixed angle of attack that is to be used by vehicle
    double fixedAngleOfAttack_;
};


//! Class containg the Pagmo-compatible formulation of the shape optimization problem.
struct ShapeOptimizationProblem
{
public:

    //! Constructor for the problem class
    /*!
     * Constructor for the problem class
     * \param bodyMap List of body objects
     * \param integratorSettings Settings for numerical integrator
     * \param propagatorSettings Settings for propagation of translational state
     * \param vehicleDensity Density of vehicle, used to compute mass from vehicle shape
     */
    ShapeOptimizationProblem(
            const simulation_setup::NamedBodyMap bodyMap,
            const std::shared_ptr< IntegratorSettings< > > integratorSettings,
            const std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings,
            const std::vector< std::pair< double, double > >& decisionVariableRange,
            const double vehicleDensity ):
        bodyMap_( bodyMap ), integratorSettings_( integratorSettings ),propagatorSettings_( propagatorSettings ),
        vehicleDensity_( vehicleDensity )
    {
        std::vector< double > boxBoundMinima, boxBoundMaxima;
        for( unsigned int i = 0; i < decisionVariableRange.size( ); i++ )
        {
            boxBoundMinima.push_back( decisionVariableRange.at( i ).first );
            boxBoundMaxima.push_back( decisionVariableRange.at( i ).second );
        }
        boxBounds_ = std::make_pair( boxBoundMinima, boxBoundMaxima );
    }

    //! Default constructor
    ShapeOptimizationProblem( ){ }

    ~ShapeOptimizationProblem( ){ }

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

    //! Function to compute propagate the dynamics of the capsule defined by the decisionVariables
    /*!
     *  Function to compute propagate the dynamics of the capsule defined by the decisionVariables. This function updates
     *  all relevant settings and properties to the new values of these parameters.
     *
     *  NOTE: Presently no fitness is computed, this must be modified during the group assignment
     *
     *  \param thrustParameters Values of parameters defining the shape and orientation of capsule (see main function)
     *  \return Fitness (undefined)
     */
    std::vector< double > fitness( const std::vector< double >& decisionVariables ) const;


    std::vector< double > getLastRunObjectives( )
    {
        return objectives_;
    }

    std::vector< double > getLastRunConstraints( )
    {
        return constraints_;
    }

    std::pair< std::vector< double >, std::vector< double > > get_bounds( ) const
    {
        return boxBounds_;
    }

private:

    void computeObjectivesAndConstraints( const std::vector< double >& decisionVariables ) const;

    //! Variable holding the body map for the simulation
    mutable simulation_setup::NamedBodyMap bodyMap_;

    //! Object holding the integrator settings
    std::shared_ptr< IntegratorSettings< > > integratorSettings_;

    //! Object holding the propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings_;

    //! Density of vehicle, used to compute mass from vehicle shape
    double vehicleDensity_;

    //! Object holding the dynamics simulator, as created during last call of fitness function
    mutable std::shared_ptr<SingleArcDynamicsSimulator< > > dynamicsSimulator_;

    //! List of objective values, as computed during last call of fitness function
    mutable std::vector< double > objectives_;

    //! List of constraint values, as computed during last call of fitness function
    mutable std::vector< double > constraints_;

    //! List of upper and lower box-bounds for decision variables.
    std::pair< std::vector< double >, std::vector< double > > boxBounds_;

};

}
}
