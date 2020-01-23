/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "shapeOptimization.h"

using namespace tudat_applications::PropagationOptimization2020;

/*!
 *  Function that creates an aerodynamic database for a capsule, based on a set of shape parameters
 *  The Capsule shape consists of four separate geometrical components: a sphere segment for the nose, a torus segment for the
 *  shoulder/edge, a conical frustum for the rear body, and a sphere segment for the rear cap (see Dirkx and Mooij, 2016).
 *  The code used in this function discretizes these surfaces into a structured mesh of quadrilateral panels. The parameters
 *  numberOfPoints and numberOfLines define the number of discretization points (for each part) in both independent directions
 *  (lengthwise and circumferential).
 */
std::shared_ptr< HypersonicLocalInclinationAnalysis > getCapsuleCoefficientInterface(
        const std::shared_ptr< geometric_shapes::Capsule > capsule,
        const std::string directory,
        const std::string filePrefix,
        const bool useNewtonianMethodForAllPanels = true )
{

    // Define settings for surface discretization of capsule
    std::vector< int > numberOfLines;
    std::vector< int > numberOfPoints;
    numberOfLines.resize( 4 );
    numberOfPoints.resize( 4 );
    numberOfLines[ 0 ] = 31;
    numberOfPoints[ 0 ] = 31;
    numberOfLines[ 1 ] = 31;
    numberOfPoints[ 1 ] = 31;
    numberOfLines[ 2 ] = 31;
    numberOfPoints[ 2 ] = 31;
    numberOfLines[ 3 ] = 11;
    numberOfPoints[ 3 ] = 11;

    // DO NOT CHANGE THESE (setting to true will turn parts of the vehicle 'inside out')
    std::vector< bool > invertOrders;
    invertOrders.resize( 4 );
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

    // Create aerodynamic database
    std::shared_ptr< HypersonicLocalInclinationAnalysis > hypersonicLocalInclinationAnalysis =
            std::make_shared< HypersonicLocalInclinationAnalysis >(
                independentVariableDataPoints, capsule, numberOfLines, numberOfPoints,
                invertOrders, selectedMethods, PI * std::pow( capsule->getMiddleRadius( ), 2.0 ),
                capsule->getMiddleRadius( ), momentReference, false );

    // Save vehicle mesh to a file
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
            const double fixedAngleOfAttack ):bodyMap_( bodyMap ), fixedAngleOfAttack_( fixedAngleOfAttack )
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

// Constructor
ShapeOptimizationProblem::ShapeOptimizationProblem()
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    outputPath_ = tudat_applications::getOutputPath( "ShapeOptimization" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            SIMULATION SETTINGS            /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

    // Vehicle properties
    double vehicleDensity = 250.0;

    // DEFINE PROBLEM INDEPENDENT VARIABLES HERE:
    std::vector< double > shapeParameters =
    { 8.148730872315355, 2.720324489288032, 0.2270385167794302, -0.4037530896422072, 0.2781438040896319, 0.4559143679738996 };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set simulation start epoch.
    simulationStartEpoch_ = 0.0;

    // Define simulation body settings.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate );
    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    bodySettings[ "Earth" ]->ephemerisSettings->resetFrameOrientation( "J2000" );

    // Create Earth object
    bodyMap_ = simulation_setup::createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create vehicle objects.
    bodyMap_[ "Capsule" ] = std::make_shared< simulation_setup::Body >( );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap_, "Earth", "J2000" );

    double limitLength = ( shapeParameters[ 1 ] - shapeParameters[ 4 ] * ( 1.0 - std::cos( shapeParameters[ 3 ] ) ) ) /
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

    bodyMap_[ "Capsule" ]->setConstantBodyMass(
                capsule->getVolume( ) * vehicleDensity );

    // Create vehicle aerodynamic coefficients
    bodyMap_[ "Capsule" ]->setAerodynamicCoefficientInterface(
                getCapsuleCoefficientInterface( capsule, outputPath_, "output_", true ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;


    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfCapsule;
    accelerationsOfCapsule[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfCapsule[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationMap[ "Capsule" ] = accelerationsOfCapsule;

    bodiesToPropagate_.push_back( "Capsule" );
    centralBodies_.push_back( "Earth" );

    // Create acceleration models
    accelerationModelMap_ = createAccelerationModelsMap(
                bodyMap_, accelerationMap, bodiesToPropagate_, centralBodies_ );

    std::shared_ptr< CapsuleAerodynamicGuidance > capsuleGuidance =
            std::make_shared< CapsuleAerodynamicGuidance >( bodyMap_, shapeParameters.at( 5 ) );
    setGuidanceAnglesFunctions( capsuleGuidance, bodyMap_.at( "Capsule" ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Convert capsule state from spherical elements to Cartesian elements.
    systemInitialState_ = convertSphericalOrbitalToCartesianState( capsuleSphericalEntryState );

    // Convert the state to the global (inertial) frame.
    systemInitialState_ = transformStateToGlobalFrame(
                systemInitialState_, simulationStartEpoch_, bodyMap_.at( "Earth" )->getRotationalEphemeris( ) );

}

void ShapeOptimizationProblem::setPropagatorIntegrator(TranslationalPropagatorType propagatorType, IntegratorType integratorType)
{
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

    // Create propagation settings
    propagatorSettings_ = std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies_, accelerationModelMap_, bodiesToPropagate_, systemInitialState_,
                terminationSettings, propagatorType, dependentVariablesToSave );


    const double initialStepSize = 1.0; // seconds
    const double minimumStepSize = 0.1; // seconds
    const double maximumStepSize = 10.0; // seconds
    const double relErrorTol = 1.0E-12;
    const double absErrorTol = 1.0E-12;
    const int    maxNumberSteps = 8;

    switch(integratorType)
    {
        case Euler:
            integratorSettings_ = std::make_shared< IntegratorSettings< > >( euler, simulationStartEpoch_, initialStepSize );
            break;
        case RK4:
            integratorSettings_ = std::make_shared< IntegratorSettings< > >( rungeKutta4, simulationStartEpoch_, initialStepSize );
            break;
        case RKF45:
            integratorSettings_ = std::make_shared< RungeKuttaVariableStepSizeSettings< > >( simulationStartEpoch_, initialStepSize, RungeKuttaCoefficients::rungeKuttaFehlberg45,
                                                                                            minimumStepSize, maximumStepSize, relErrorTol, absErrorTol );
            break;
        case RKF56:
            integratorSettings_ = std::make_shared< RungeKuttaVariableStepSizeSettings< > >( simulationStartEpoch_, initialStepSize, RungeKuttaCoefficients::rungeKuttaFehlberg56,
                                                                                            minimumStepSize, maximumStepSize, relErrorTol, absErrorTol );
            break;
        case RKF78:
            integratorSettings_ = std::make_shared< RungeKuttaVariableStepSizeSettings< > >( simulationStartEpoch_, initialStepSize, RungeKuttaCoefficients::rungeKuttaFehlberg78,
                                                                                        minimumStepSize, maximumStepSize, relErrorTol, absErrorTol );
            break;
        case RK87:
            integratorSettings_ = std::make_shared< RungeKuttaVariableStepSizeSettings< > >( simulationStartEpoch_, initialStepSize, RungeKuttaCoefficients::rungeKutta87DormandPrince,
                                                                                    minimumStepSize, maximumStepSize, relErrorTol, absErrorTol );
            break;
        case ABM:
            // ABM
            integratorSettings_ = std::make_shared< AdamsBashforthMoultonSettings< > >( simulationStartEpoch_, initialStepSize, minimumStepSize, maximumStepSize, relErrorTol,
                                                                                       absErrorTol, 6, 11);
            break;
        case BS:
            // BS
            integratorSettings_ = std::make_shared< BulirschStoerIntegratorSettings< > >( simulationStartEpoch_, initialStepSize, ExtrapolationMethodStepSequences::bulirsch_stoer_sequence,
                                                                                         maxNumberSteps, minimumStepSize, maximumStepSize, relErrorTol, absErrorTol);
            break;
    }

}


// The const identifier is really annoying (but required by Pagmo), as we cannot set class member fields or call non-const functions to set them for us
std::vector< double > ShapeOptimizationProblem::fitness( std::vector< double >& x) const
{

    // Use parameter x in some way

    // Create simulation object and propagate dynamics
    SingleArcDynamicsSimulator< > dynamicsSimulator{ bodyMap_, integratorSettings_, propagatorSettings_ };

    // Call separate class methods to return the state and dependent variable maps
    returnPropagatedStateHistory( dynamicsSimulator.getEquationsOfMotionNumericalSolution( ) );
    returnDependentVariableHistory( dynamicsSimulator.getDependentVariableHistory( ) );

    // Return fitness; for now just return {0.0}
    return {0.0};

}
