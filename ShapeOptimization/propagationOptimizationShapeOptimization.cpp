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

//! Function that creates an aerodynamic database for a capsule, based on a set of shape parameters
//! (see main function documentation)
std::shared_ptr< HypersonicLocalInclinationAnalysis > getCapsuleCoefficientInterface(
        const std::shared_ptr< geometric_shapes::Capsule > capsule,
        const std::string directory,
        const std::string filePrefix,
        const bool useNewtonianMethodForAllPanels = true )
{

    // Define geometry settings
    std::vector< int > numberOfLines;
    std::vector< int > numberOfPoints;
    std::vector< bool > invertOrders;
    numberOfLines.resize( 4 );
    numberOfPoints.resize( 4 );
    invertOrders.resize( 4 );
    numberOfLines[ 0 ] = 31;
    numberOfPoints[ 0 ] = 31;
    numberOfLines[ 1 ] = 31;
    numberOfPoints[ 1 ] = 31;
    numberOfLines[ 2 ] = 31;
    numberOfPoints[ 2 ] = 31;
    numberOfLines[ 3 ] = 11;
    numberOfPoints[ 3 ] = 11;
    invertOrders[ 0 ] = 0;
    invertOrders[ 1 ] = 0;
    invertOrders[ 2 ] = 0;
    invertOrders[ 3 ] = 0;

    // Define moment reference point
    Eigen::Vector3d momentReference;
    momentReference( 0 ) = -0.6624;
    momentReference( 1 ) = 0.0;
    momentReference( 2 ) = 0.1369;

    // Define independent variable values
    std::vector< std::vector< double > > independentVariableDataPoints;
    independentVariableDataPoints.resize( 3 );
    independentVariableDataPoints[ 0 ] = getDefaultHypersonicLocalInclinationMachPoints( "Full" );
    std::vector< double > angleOfAttackPoints;
    angleOfAttackPoints.resize( 15 );
    for ( int i = 0; i < 15; i++ )
    {
        angleOfAttackPoints[ i ] = static_cast< double >( i - 6 ) * 5.0 * PI / 180.0;
    }
    independentVariableDataPoints[ 1 ] = angleOfAttackPoints;
    independentVariableDataPoints[ 2 ] =
            getDefaultHypersonicLocalInclinationAngleOfSideslipPoints( );

    // Define local inclination methods to use
    std::vector< std::vector< int > > selectedMethods;
    selectedMethods.resize( 2 );
    selectedMethods[ 0 ].resize( 4 );
    selectedMethods[ 1 ].resize( 4 );
    if( !useNewtonianMethodForAllPanels )
    {
        selectedMethods[ 0 ][ 0 ] = 1;
        selectedMethods[ 0 ][ 1 ] = 5;
        selectedMethods[ 0 ][ 2 ] = 5;
        selectedMethods[ 0 ][ 3 ] = 1;
        selectedMethods[ 1 ][ 0 ] = 6;
        selectedMethods[ 1 ][ 1 ] = 3;
        selectedMethods[ 1 ][ 2 ] = 3;
        selectedMethods[ 1 ][ 3 ] = 3;
    }
    else
    {
        selectedMethods[ 0 ][ 0 ] = 0;
        selectedMethods[ 0 ][ 1 ] = 0;
        selectedMethods[ 0 ][ 2 ] = 0;
        selectedMethods[ 0 ][ 3 ] = 0;
        selectedMethods[ 1 ][ 0 ] = 0;
        selectedMethods[ 1 ][ 1 ] = 0;
        selectedMethods[ 1 ][ 2 ] = 0;
        selectedMethods[ 1 ][ 3 ] = 0;
    }

    std::shared_ptr< HypersonicLocalInclinationAnalysis > hypersonicLocalInclinationAnalysis =
            std::make_shared< HypersonicLocalInclinationAnalysis >(
                independentVariableDataPoints, capsule, numberOfLines, numberOfPoints,
                invertOrders, selectedMethods, PI * std::pow( capsule->getMiddleRadius( ), 2.0 ),
                capsule->getMiddleRadius( ), momentReference, false );

    aerodynamics::saveVehicleMeshToFile(
                hypersonicLocalInclinationAnalysis, directory, filePrefix );

    // Create analysis object and capsule database.
    return  hypersonicLocalInclinationAnalysis;
}

namespace tudat
{

//! Class to set the aerodynamic angles of the capsule (default: all angles 0)
class CapsuleAerodynamicGuidance: public aerodynamics::AerodynamicGuidance
{
public:

    //! Constructor
    CapsuleAerodynamicGuidance(
            const NamedBodyMap bodyMap,
            const double fixedAngleOfAttack ):bodyMap_( bodyMap )
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

}

/*!
 *   This function computes the entry trajectory of a capsule, where the shape of the capsule is used to determine the vehicle's
 *   aerodynamic force and moment coefficients. The aerodynamic coefficients are based on local inclination methods, and computed
 *   by an object of the HypersonicLocalInclinationAnalysis class.
 *
 *   The present code uses only an aerodynamic force, and a point-mass Earth gravitational acceleration. A class
 *   CapsuleAerodynamicGuidance has been provided, which is currently has no direct functionality: it sets the aerodynamic angles
 *   (attack, sideslip, bank) to 0 degrees. This can (and should) be overridden by the user in favor of something more realistic
 *
 *   Key outputs:
 *
 *   propagatedStateHistory Numerically propagated Cartesian state
 *   dependentVariableHistory Dependent variables (default none) saved during the state propagation of the entry capsule
 *
 *   Input parameters:
 *
 *   shapeParameters: A vector defining the shape of the capsule as follows: Nose radius, Middle radius, Rear length, Rear angle
 *      Side radius, Constant Angle of Attack (see Dirkx and Mooij, 2018 for more details)
 *
 */
//! Execute propagation of orbits of Capsule during entry.
int main( )
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );


    std::string outputPath = tudat_applications::getOutputPath( "ShapeOptimization" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            SIMULATION SETTINGS            /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    double vehicleDensity = 250.0;

    // Set spherical elements for Capsule.
    Eigen::Vector6d capsuleSphericalEntryState;
    capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
            spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
    capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) =
            unit_conversions::convertDegreesToRadians( 0.0 );
    capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) =
            unit_conversions::convertDegreesToRadians( 68.75 );
    capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 7.83E3;
    capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
            unit_conversions::convertDegreesToRadians( -1.5 );
    capsuleSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) =
            unit_conversions::convertDegreesToRadians( 34.37 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 1.0;

    // Define simulation body settings.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate );
    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    bodySettings[ "Earth" ]->ephemerisSettings->resetFrameOrientation( "J2000" );

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create vehicle objects.
    bodyMap[ "Capsule" ] = std::make_shared< simulation_setup::Body >( );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "Earth", "J2000" );

    // DEFINE PROBLEM INDEPENDENT VARIABLES HERE:
    std::vector< double > shapeParameters;

    double limitLength =
            ( shapeParameters[ 1 ] - shapeParameters[ 4 ] * ( 1.0 - std::cos( shapeParameters[ 3 ] ) ) ) /
            std::tan( -shapeParameters[ 3 ] );
    if( shapeParameters[ 2 ] >= limitLength - 0.01 )
    {
        shapeParameters[ 2 ] = limitLength -0.01;
    }

    // Create capsule.
    std::shared_ptr< geometric_shapes::Capsule > capsule
            = std::make_shared< geometric_shapes::Capsule >(
                shapeParameters[ 0 ], shapeParameters[ 1 ], shapeParameters[ 2 ],
            shapeParameters[ 3 ], shapeParameters[ 4 ] );

    bodyMap[ "Capsule" ]->setConstantBodyMass(
                capsule->getVolume( ) * vehicleDensity );

    // Create vehicle aerodynamic coefficients
    bodyMap[ "Capsule" ]->setAerodynamicCoefficientInterface(
                getCapsuleCoefficientInterface( capsule, outputPath, "output_", true ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfCapsule;
    accelerationsOfCapsule[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfCapsule[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationMap[ "Capsule" ] = accelerationsOfCapsule;

    bodiesToPropagate.push_back( "Capsule" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    std::shared_ptr< CapsuleAerodynamicGuidance > capsuleGuidance =
            std::make_shared< CapsuleAerodynamicGuidance >( bodyMap, shapeParameters.at( 5 ) );
    setGuidanceAnglesFunctions( capsuleGuidance, bodyMap.at( "Capsule" ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Convert capsule state from spherical elements to Cartesian elements.
    Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState(
                capsuleSphericalEntryState );

    // Convert the state to the global (inertial) frame.
    systemInitialState = transformStateToGlobalFrame(
                systemInitialState, simulationStartEpoch, bodyMap.at( "Earth" )->getRotationalEphemeris( ) );

    // Define termination conditions
    std::vector< std::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;

    std::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
            std::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "Capsule", "Earth" );
    terminationSettingsList.push_back(
                std::make_shared< PropagationDependentVariableTerminationSettings >(
                    terminationDependentVariable, 25.0E3, true ) );
    terminationSettingsList.push_back( std::make_shared< PropagationTimeTerminationSettings >(
                                           24.0 * 3600.0 ) );
    std::shared_ptr< PropagationTerminationSettings > terminationSettings = std::make_shared<
            PropagationHybridTerminationSettings >( terminationSettingsList, true );

    // Define dependent variables
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          altitude_dependent_variable, "Capsule", "Earth" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          airspeed_dependent_variable, "Capsule", "Earth" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          aerodynamic_force_coefficients_dependent_variable, "Capsule" ) );

    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

    // Create propagation settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                terminationSettings, cowell, dependentVariablesToSave );
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >( rungeKutta4, simulationStartEpoch, fixedStepSize );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings );

    std::map< double, Eigen::VectorXd > propagatedStateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );


    input_output::writeDataMapToTextFile( propagatedStateHistory, "stateHistory.dat", outputPath );
    input_output::writeDataMapToTextFile( dependentVariableHistory, "dependentVariables.dat", outputPath );

    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
