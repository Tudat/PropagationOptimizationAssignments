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

Eigen::Vector6d getInitialState( double simulationStartEpoch, simulation_setup::NamedBodyMap bodyMap )
{
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

    // Set initial inertial Cartesian state
    Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState( capsuleSphericalEntryState );
    systemInitialState = transformStateToGlobalFrame(
                systemInitialState, simulationStartEpoch, bodyMap.at( "Earth" )->getRotationalEphemeris( ) );

    return systemInitialState;

}

std::shared_ptr< PropagationTerminationSettings > getPropagationTerminationSettings()
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
    return std::make_shared< PropagationHybridTerminationSettings >( terminationSettingsList, true );
}

std::shared_ptr< DependentVariableSaveSettings > getDependentVariableSaveSettings()
{
    // Define dependent variables
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          mach_number_dependent_variable, "Capsule" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          altitude_dependent_variable, "Capsule", "Earth" ) );
    return std::make_shared< DependentVariableSaveSettings >( dependentVariablesList, false );
}

std::vector< std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > > generateBenchmarks(
        double simulationStartEpoch, simulation_setup::NamedBodyMap bodyMap,
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > benchmarkPropagatorSettings,
        std::vector< double > shapeParameters, std::string outputPath )
{
    // Create integrator settings
    std::shared_ptr< IntegratorSettings< > > benchmarkIntegratorSettings;
    benchmarkIntegratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                simulationStartEpoch, 0.2, RungeKuttaCoefficients::rungeKutta87DormandPrince,
                0.2, 0.2,
                std::numeric_limits< double >::infinity( ), std::numeric_limits< double >::infinity( ) );

    ShapeOptimizationProblem probBenchmarkFirst{ bodyMap, benchmarkIntegratorSettings, benchmarkPropagatorSettings };
    probBenchmarkFirst.fitness( shapeParameters );


    benchmarkIntegratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                simulationStartEpoch, 0.4, RungeKuttaCoefficients::rungeKutta87DormandPrince,
                0.4, 0.4,
                std::numeric_limits< double >::infinity( ), std::numeric_limits< double >::infinity( ) );
    ShapeOptimizationProblem probBenchmarkSecond{ bodyMap, benchmarkIntegratorSettings, benchmarkPropagatorSettings };
    probBenchmarkSecond.fitness( shapeParameters );

    std::map< double, Eigen::VectorXd > firstBenchmarkResults = probBenchmarkFirst.getLastRunPropagatedStateHistory( );
    std::map< double, Eigen::VectorXd > secondBenchmarkResults = probBenchmarkSecond.getLastRunPropagatedStateHistory( );

    input_output::writeDataMapToTextFile( firstBenchmarkResults,
                                          "benchmark1.dat", outputPath );
    input_output::writeDataMapToTextFile( secondBenchmarkResults,
                                          "benchmark2.dat", outputPath );

    // Create 8th-order Lagrange interpolator
    std::shared_ptr< InterpolatorSettings > interpolatorSettings = std::make_shared< LagrangeInterpolatorSettings >( 8 );
    std::vector< std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > > interpolators;

    interpolators.push_back( createOneDimensionalInterpolator( firstBenchmarkResults, interpolatorSettings ) );
    interpolators.push_back( createOneDimensionalInterpolator( secondBenchmarkResults, interpolatorSettings ) );

    return interpolators;
}

int main()
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );
    std::string outputPath = tudat_applications::getOutputPath( "ShapeOptimization" );

    std::vector< double > shapeParameters =
    { 8.148730872315355, 2.720324489288032, 0.2270385167794302, -0.4037530896422072, 0.2781438040896319, 0.4559143679738996 };


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
    addCapsuleToBodyMap( bodyMap, shapeParameters );

    // Create initial Cartesian state from spherical elements
    Eigen::Vector6d systemInitialState = getInitialState( simulationStartEpoch, bodyMap );

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

    //    // Create acceleration models
    //    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
    //                bodyMap, accelerationSettingsMap, bodiesToPropagate, centralBodies );
    //    std::shared_ptr< CapsuleAerodynamicGuidance > capsuleGuidance =
    //            std::make_shared< CapsuleAerodynamicGuidance >( bodyMap, shapeParameters.at( 5 ) );
    //    setGuidanceAnglesFunctions( capsuleGuidance, bodyMap.at( "Capsule" ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::shared_ptr< PropagationTerminationSettings > terminationSettings = getPropagationTerminationSettings();
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave = getDependentVariableSaveSettings();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////  IF DESIRED, GENERATE BENCHMARK                            ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool generateAndCompareToBenchmark = true;
    if( generateAndCompareToBenchmark )
    {
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > benchmarkPropagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies, accelerationSettingsMap, bodiesToPropagate, systemInitialState,
                    terminationSettings, cowell, dependentVariablesToSave );
        std::vector< std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > > benchmarkInterpolators =
                generateBenchmarks(simulationStartEpoch, bodyMap, benchmarkPropagatorSettings, shapeParameters, outputPath);

    }

    // TO BE DONE: use benchmark compute difference w.r.t. all other runs

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             RUN SIMULATION FOR VARIOUS SETTINGS            ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define list of propagators (for convenience)
    std::vector< TranslationalPropagatorType > propagatorTypes =
    { cowell, encke, gauss_keplerian, gauss_modified_equinoctial,
      unified_state_model_quaternions, unified_state_model_modified_rodrigues_parameters,
      unified_state_model_exponential_map };

    // Define list of multi-stage integrators (for convenience)
    std::vector< RungeKuttaCoefficients::CoefficientSets > multiStageTypes =
    { RungeKuttaCoefficients::rungeKuttaFehlberg45,
      RungeKuttaCoefficients::rungeKuttaFehlberg56,
      RungeKuttaCoefficients::rungeKuttaFehlberg78,
      RungeKuttaCoefficients::rungeKutta87DormandPrince };

    // TODO: map propagators and integrators to numbers and output to file

    // Define number of settings to use
    int numberOfPropagators = 7;
    int numberOfIntegrators = 5;
    int numberOfIntegratorStepSizeSettings = 4;
    for( int i = 0; i < numberOfPropagators; i++ )
    {
        // Create propagator settings
        TranslationalPropagatorType propagatorType = propagatorTypes.at( i );
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies, accelerationSettingsMap, bodiesToPropagate, systemInitialState,
                    terminationSettings, propagatorType, dependentVariablesToSave );
        //        propagatorSettings->resetIntegratedStateModels( bodyMap );

        // Iterate over all types of integrators
        for( int j = 0; j < numberOfIntegrators; j++ )
        {
            // Change number of integrator settings for RK4
            if( j >= 4 )
            {
                numberOfIntegratorStepSizeSettings = 6;
            }

            // Iterate over all tolerances/step sizes
            for( int k = 0; k < numberOfIntegratorStepSizeSettings; k++ )
            {
                // Print status
                std::cout<<"Current run "<<i<<" "<<j<<" "<<k<<std::endl;

                // Define integrator settings
                std::shared_ptr< IntegratorSettings< > > integratorSettings;

                // Create integrator settings (multi-stage variable-step)
                if( j < 4 )
                {
                    // Extract integrator type and tolerance for current run
                    RungeKuttaCoefficients::CoefficientSets currentCoefficientSet = multiStageTypes.at( j );
                    double currentTolerance = std::pow( 10.0, -14 + k );

                    // CONVERT TO FUNCTION (j, k) -> INTEGRATOR SETTINGS OBJECT

                    // Create integrator settings
                    integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                                simulationStartEpoch, 1.0, currentCoefficientSet,
                                std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                                currentTolerance, currentTolerance );

                }
                // Create integrator settings (multi-stage fixed-step)
                else
                {
                    // Create integrator settings
                    double timeStep = 0.1 * std::pow( 2, k );
                    integratorSettings =
                            std::make_shared< IntegratorSettings< > >
                            ( rungeKutta4, simulationStartEpoch, timeStep );
                }


                // Construct problem
                ShapeOptimizationProblem prob{ bodyMap, integratorSettings, propagatorSettings };

                // Propagate trajectory using defined settings
                prob.fitness( shapeParameters );

                // Save state and dependent variable results to file
                input_output::writeDataMapToTextFile( prob.getLastRunPropagatedStateHistory( ),
                                                      "stateHistory_test_"
                                                      + std::to_string( i ) + "_"
                                                      + std::to_string( j ) + "_"
                                                      + std::to_string( k ) + "_" + ".dat", outputPath );
                input_output::writeDataMapToTextFile( prob.getLastRunDependentVariableHistory( ),
                                                      "dependentVariables_test_"
                                                      + std::to_string( i ) + "_"
                                                      + std::to_string( j ) + "_"
                                                      + std::to_string( k ) + "_" + ".dat", outputPath );
            }
        }
    }

}
