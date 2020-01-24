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

int main()
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );
    std::string outputPath = tudat_applications::getOutputPath( "ShapeOptimization" );

    std::vector< double > shapeParameters =
    { 8.148730872315355, 2.720324489288032, 0.2270385167794302, -0.4037530896422072, 0.2781438040896319, 0.4559143679738996 };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            BASIC SIMULATION SETTINGS            ///////////////////////////////////////////////
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


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set simulation start epoch.
    double simulationStartEpoch = 0.0;

    // Define simulation body settings.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate );
    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    bodySettings[ "Earth" ]->ephemerisSettings->resetFrameOrientation( "J2000" );

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

    // Create vehicle objects.
    bodyMap[ "Capsule" ] = std::make_shared< simulation_setup::Body >( );
    setVehicleShapeParameters( shapeParameters, bodyMap );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "Earth", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationSettingsMap;

    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfCapsule;
    accelerationsOfCapsule[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfCapsule[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationSettingsMap[ "Capsule" ] = accelerationsOfCapsule;

    std::vector< std::string > centralBodies;
    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back( "Capsule" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationSettingsMap, bodiesToPropagate, centralBodies );
    std::shared_ptr< CapsuleAerodynamicGuidance > capsuleGuidance =
            std::make_shared< CapsuleAerodynamicGuidance >( bodyMap, shapeParameters.at( 5 ) );
    setGuidanceAnglesFunctions( capsuleGuidance, bodyMap.at( "Capsule" ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set initial inertial Cartesian state
    Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState( capsuleSphericalEntryState );
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
                                          mach_number_dependent_variable, "Capsule" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          altitude_dependent_variable, "Capsule", "Earth" ) );
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList, false );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             RUN SIMULATION FOR VARIOUS SETTINGS            ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    std::vector< TranslationalPropagatorType > propagatorTypes =
    { cowell, encke, gauss_keplerian, gauss_modified_equinoctial,
      unified_state_model_quaternions, unified_state_model_modified_rodrigues_parameters,
      unified_state_model_exponential_map };


    std::vector< RungeKuttaCoefficients::CoefficientSets > multiStageTypes =
    { RungeKuttaCoefficients::rungeKuttaFehlberg45,
            RungeKuttaCoefficients::rungeKuttaFehlberg56,
            RungeKuttaCoefficients::rungeKuttaFehlberg78,
            RungeKuttaCoefficients::rungeKutta87DormandPrince };

    int numberOfPropagators = 7;
    int numberOfIntegrators = 5;
    int numberOfIntegratorStepSizeSettings = 4;
    for( int i = 0; i < numberOfPropagators; i++ )
    {
        TranslationalPropagatorType propagatorType = propagatorTypes.at( i );
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                    terminationSettings, propagatorType, dependentVariablesToSave );

        for( int j = 0; j < numberOfIntegrators; j++ )
        {
            if( j < 4 )
            {
                numberOfIntegratorStepSizeSettings = 4;
            }
            else
            {
                numberOfIntegratorStepSizeSettings = 6;
            }
            for( int k = 0; k < numberOfIntegratorStepSizeSettings; k++ )
            {
                std::cout<<"Current run "<<i<<" "<<j<<" "<<k<<std::endl;
                std::shared_ptr< IntegratorSettings< > > integratorSettings;
                if( j < 4 )
                {
                    integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                                        simulationStartEpoch, 1.0,
                                        multiStageTypes.at( j ),
                                        std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                                        std::pow( 10.0, -14 + k ) );

                }
                else
                {
                    integratorSettings =
                            std::make_shared< IntegratorSettings< > >
                            ( rungeKutta4, 0.0, 0.1 * std::pow( 2, k ) );
                }


                // Construct problem
                ShapeOptimizationProblem prob{ bodyMap, integratorSettings, propagatorSettings };

                // Propagate trajectory using defined settings
                prob.fitness( shapeParameters );

                // TODO: retrieve state & dependent variable history from fitness function (const function!)
                input_output::writeDataMapToTextFile( prob.getLastRunPropagatedStateHistory( ),
                                                      "stateHistory_"
                                                      + std::to_string(i) + "_"
                                                       + std::to_string(j) + "_"
                                                       + std::to_string(k) + "_" + ".dat", outputPath );
                input_output::writeDataMapToTextFile( prob.getLastRunDependentVariableHistory( ),
                                                      "dependentVariables_"
                                                      + std::to_string(i) + "_"
                                                       + std::to_string(j) + "_"
                                                       + std::to_string(k) + "_" + ".dat", outputPath );
            }
        }
    }

}
