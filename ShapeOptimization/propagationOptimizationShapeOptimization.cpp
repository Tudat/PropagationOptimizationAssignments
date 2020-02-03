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

//! Function to retrieve the initial Cartesian state of the vehicle.
/*!
 * Function to retrieve the initial Cartesian state of the vehicle. The spherical orbital parameters are
 * first converted to Cartesian coordinates and subsequently transformed to the global frame of reference.
 * \param simulationStartEpoch The start time of the simulation in seconds.
 * \param bodyMap NamedBodyMap containing the bodies in the simulation.
 * \return Eigen Vector6d containing the system's initial state in Cartesian coordinates.
 */
Eigen::Vector6d getInitialState( double simulationStartEpoch, simulation_setup::NamedBodyMap bodyMap )
{
    // Set spherical elements for Capsule
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

    // Set initial inertial Cartesian state and convert to global frame of reference
    Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState( capsuleSphericalEntryState );
    systemInitialState = transformStateToGlobalFrame(
                systemInitialState, simulationStartEpoch, bodyMap.at( "Earth" )->getRotationalEphemeris( ) );

    return systemInitialState;

}

//! Get the propagation termination settings for the shape optimization.
/*!
 * This function returns a shared pointer to a PropagationTerminationSettings object, containing
 * altitude and time termination settings.
 * \return Shared pointer to the PropagationTerminationSettings object.
 */
std::shared_ptr< PropagationTerminationSettings > getPropagationTerminationSettings()
{
    // Define termination conditions
    std::vector< std::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;
    std::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
            std::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "Capsule", "Earth" );
    // Altitude termination condition
    terminationSettingsList.push_back(
                std::make_shared< PropagationDependentVariableTerminationSettings >(
                    terminationDependentVariable, 25.0E3, true ) );
    // Time termination condition
    terminationSettingsList.push_back( std::make_shared< PropagationTimeTerminationSettings >(
                                           24.0 * 3600.0 ) );
    return std::make_shared< PropagationHybridTerminationSettings >( terminationSettingsList, true );
}

//! Function to retrieve the integrator settings for the current run.
/*!
 * This function returns a shared pointer to an IntegratorSettings object, based on the indices passed to
 * the function. The first index is used to determine which type is selected from the vector of integrator
 * types which is passed as the last parameter. The seconds index determines what tolerance is used for the
 * variable step size integrators.
 * \param j Index specifying which kind of integrator is used: if smaller than 4, use a multi-stage,
 * variable step size RKF or RK Dormand-Prince; otherwise, use an RK4 with fixed step size.
 * \param k Index that is used to specify different tolerances for the same integrator. k = 0
 * corresponds to 10^-14, k = 1 to 10^-13, etc.
 * \param simulationStartEpoch The start time of the simulation in seconds.
 * \param multiStageTypes std::vector with Runge-Kutta coefficient sets to use.
 * \return Shared pointer to the IntegratorSettings object.
 */
std::shared_ptr< IntegratorSettings< > > getIntegratorSettings( unsigned int j, int k, double simulationStartEpoch,
                                                                std::vector< RungeKuttaCoefficients::CoefficientSets > multiStageTypes )
{
    // Create integrator settings (multi-stage variable-step)
    if( j < 4 )
    {
        // Extract integrator type and tolerance for current run
        RungeKuttaCoefficients::CoefficientSets currentCoefficientSet = multiStageTypes.at( j );
        double currentTolerance = std::pow( 10.0, -14 + k );

        // Create integrator settings
        return std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                    simulationStartEpoch, 1.0, currentCoefficientSet,
                    std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                    currentTolerance, currentTolerance );

    }
    // Create integrator settings (multi-stage fixed-step)
    else
    {
        // Create integrator settings
        double timeStep = 0.1 * std::pow( 2, k );
        return std::make_shared< IntegratorSettings< > >( rungeKutta4, simulationStartEpoch, timeStep );
    }
}

//! Function to retrieve the dependent variable save settings for the current simulation.
/*!
 * This function returns a shared pointer to a DependentVariableSaveSettings object, containing the save settings
 * to save the Mach number and altitude of the capsule.
 * \return Shared pointer to a DependentVariableSaveSettings object.
 */
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

//! Function to generate to accurate benchmarks.
/*!
 * This function runs two propagations with two different integrator settings that serve as benchmarks for
 * the nominal runs. To be able to compare these, the function returns the two interpolators pertaining
 * to the state and dependent variables of one of the benchmarks. The states are written to files, as well
 * as the difference in state and dependent variables between the two benchmarks.
 * \param simulationStartEpoch The start time of the simulation in seconds.
 * \param bodyMap NamedBodyMap containing the bodies in the simulation.
 * \param benchmarkPropagatorSettings Shared pointer to a translational propagator settings object,
 * which is used to run the benchmark propagations.
 * \param shapeParameters The vector of doubles that represents the shape parameters for the capsule.
 * \param outputPath String containing the path to the output directory.
 * \return
 */
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

    std::cout << "Running first benchmark..." << std::endl;
    probBenchmarkFirst.fitness( shapeParameters );


    benchmarkIntegratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                simulationStartEpoch, 0.4, RungeKuttaCoefficients::rungeKutta87DormandPrince,
                0.4, 0.4,
                std::numeric_limits< double >::infinity( ), std::numeric_limits< double >::infinity( ) );
    ShapeOptimizationProblem probBenchmarkSecond{ bodyMap, benchmarkIntegratorSettings, benchmarkPropagatorSettings };

    std::cout << "Running second benchmark..." << std::endl;
    probBenchmarkSecond.fitness( shapeParameters );

    std::map< double, Eigen::VectorXd > firstBenchmarkStates = probBenchmarkFirst.getLastRunPropagatedStateHistory( );
    std::map< double, Eigen::VectorXd > secondBenchmarkStates = probBenchmarkSecond.getLastRunPropagatedStateHistory( );

    std::map< double, Eigen::VectorXd > firstBenchmarkDependent = probBenchmarkFirst.getLastRunDependentVariableHistory( );
    std::map< double, Eigen::VectorXd > secondBenchmarkDependent = probBenchmarkSecond.getLastRunDependentVariableHistory( );
\
    // Put the benchmark data in a separate directory
    outputPath.append("/benchmarks/");

    // Write the state maps of both benchmarks to files
    input_output::writeDataMapToTextFile( firstBenchmarkStates,
                                          "benchmark1.dat", outputPath );
    input_output::writeDataMapToTextFile( secondBenchmarkStates,
                                          "benchmark2.dat", outputPath );

    // Create 8th-order Lagrange interpolator
    std::shared_ptr< InterpolatorSettings > interpolatorSettings = std::make_shared< LagrangeInterpolatorSettings >( 8 );
    std::vector< std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > > interpolators;

    std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > firstStatesInterpolator =
            createOneDimensionalInterpolator( firstBenchmarkStates, interpolatorSettings );
    std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > secondStatesInterpolator =
            createOneDimensionalInterpolator( secondBenchmarkStates, interpolatorSettings );
    std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > firstDependentInterpolator =
            createOneDimensionalInterpolator( firstBenchmarkDependent, interpolatorSettings );
    std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > secondDependentInterpolator =
            createOneDimensionalInterpolator( secondBenchmarkDependent, interpolatorSettings );

    std::cout << "Calculating the difference between the benchmarks..." << std::endl;
    std::map< double, Eigen::VectorXd > statesDifference;
    std::map< double, Eigen::VectorXd > dependentVariablesDifference;

    // Calculate difference between two benchmarks, both for the states and the dependent variables
    for( auto iterator = firstBenchmarkStates.begin(); iterator != firstBenchmarkStates.end(); iterator++ )
    {
        statesDifference[ iterator->first ] = secondStatesInterpolator->interpolate( iterator->first ) - iterator->second;
    }
    for( auto iterator = firstBenchmarkDependent.begin(); iterator != firstBenchmarkDependent.end(); iterator++ )
    {
        dependentVariablesDifference[ iterator->first ] = secondDependentInterpolator->interpolate( iterator->first ) - iterator->second;
    }

    // Write the difference in state and dependent variables between the two benchmarks to files
    input_output::writeDataMapToTextFile( statesDifference,
                                          "benchmarkStateDifference.dat", outputPath );
    input_output::writeDataMapToTextFile( dependentVariablesDifference,
                                          "benchmarkDependentDifference.dat", outputPath );

    // Return the interpolators for the first benchmark (can be changed to the second if needed)
    interpolators.push_back( firstStatesInterpolator );
    interpolators.push_back( firstDependentInterpolator );

    return interpolators;
}

//! Function to create two files, specifying which index belongs to which integrator and propagator.
/*!
 * \param outputPath String containing the path to the output directory.
 */
void writePropagatorIntegratorIndicesToFile( std::string outputPath )
{
    std::map< unsigned int, std::string > propagatorIndices;
    std::map< unsigned int, std::string > integratorIndices;

    propagatorIndices[0] = "Cowell";
    propagatorIndices[1] = "Encke";
    propagatorIndices[2] = "Gauss Keplerian";
    propagatorIndices[3] = "Gauss Modified Equinoctial";
    propagatorIndices[4] = "USM Quaternions";
    propagatorIndices[5] = "USM MRP";
    propagatorIndices[6] = "USM Exponential Map";

    // TODO: find another way to make this more robust/dynamic
    integratorIndices[0] = "RKF45";
    integratorIndices[1] = "RKF56";
    integratorIndices[2] = "RKF78";
    integratorIndices[3] = "RK87 DP";
    integratorIndices[4] = "RK4";

    input_output::writeDataMapToTextFile( propagatorIndices, "propagator_indices.txt", outputPath );
    input_output::writeDataMapToTextFile( integratorIndices, "integrator_indices.txt", outputPath );
}

/*!
 *   This function computes the dynamics of a capsule re-entering the atmosphere of the Earth, starting 120 km above the surface
 *   of the planet, with a velocity of 7.83 km/s. This propagation assumes only point mass gravity and aerodynamic acceleration, providing
 *   a baseline for further investigation. The shape of the capsule can be tweaked by the shape parameters to find the optimum
 *   shape of the vehicle.
 *
 *   The trajectory of the capsule is heavily dependent on the shape of the vehicle, which is determined by the five shape parameters
 *   that form the input for the simulation (see below). These parameters are used to compute the aerodynamic accelerations on the
 *   vehicle (see also [PUBLICATION BY DOMINIC DIRKX, TODO]).
 *
 *   The propagation is terminated as soon as one of the following conditions is met:
 *
 *   - Altitude < 25 km
 *   - Propagation time > 24 hr
 *
 *   Key outputs:
 *
 *   propagatedStateHistory Numerically propagated Cartesian state
 *   dependentVariableHistory Dependent variables saved during the state propagation of the ascent *
 *
 *   Input parameters:
 *
 *   shapeParameters: Vector contains the following:
 *
 *   - Entry 0: TO BE SPECIFIED
 *   - Entry 1: TO BE SPECIFIED
 *   - Entry 2: TO BE SPECIFIED
 *   - Entry 3: TO BE SPECIFIED
 *   - Entry 4: TO BE SPECIFIED
 *   - Entry 5: TO BE SPECIFIED
 */
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
    std::vector< std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > > benchmarkInterpolators;
    if( generateAndCompareToBenchmark )
    {
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > benchmarkPropagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies, accelerationSettingsMap, bodiesToPropagate, systemInitialState,
                    terminationSettings, cowell, dependentVariablesToSave );
         benchmarkInterpolators = generateBenchmarks(simulationStartEpoch, bodyMap, benchmarkPropagatorSettings,
                                                     shapeParameters, outputPath);

    }

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

    writePropagatorIntegratorIndicesToFile( outputPath );

    // Define number of settings to use
    const unsigned int numberOfPropagators = 7;
    const unsigned int numberOfIntegrators = 5;
    int numberOfIntegratorStepSizeSettings = 4;
    for( unsigned int i = 0; i < numberOfPropagators; i++ )
    {
        // Create propagator settings
        TranslationalPropagatorType propagatorType = propagatorTypes.at( i );
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies, accelerationSettingsMap, bodiesToPropagate, systemInitialState,
                    terminationSettings, propagatorType, dependentVariablesToSave );
        //        propagatorSettings->resetIntegratedStateModels( bodyMap );

        // Iterate over all types of integrators
        for( unsigned int j = 0; j < numberOfIntegrators; j++ )
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

                outputPath = tudat_applications::getOutputPath( "ShapeOptimization/" + std::to_string( i ) + "/" +
                                                                std::to_string( j ) + "/" + std::to_string( k )
                                                                + "/" );

                // Define integrator settings
                std::shared_ptr< IntegratorSettings< > > integratorSettings =
                        getIntegratorSettings( j, k, simulationStartEpoch, multiStageTypes );

                // Construct problem
                ShapeOptimizationProblem prob{ bodyMap, integratorSettings, propagatorSettings };

                // Propagate trajectory using defined settings
                prob.fitness( shapeParameters );

                std::map< double, Eigen::VectorXd> stateHistory = prob.getLastRunPropagatedStateHistory( );
                std::map< double, Eigen::VectorXd> dependentVariableHistory = prob.getLastRunDependentVariableHistory( );

                // Save state and dependent variable results to file
                input_output::writeDataMapToTextFile( stateHistory,
                                                      "stateHistory_test_"
                                                      + std::to_string( i ) + "_"
                                                      + std::to_string( j ) + "_"
                                                      + std::to_string( k ) + ".dat", outputPath );
                input_output::writeDataMapToTextFile( dependentVariableHistory,
                                                      "dependentVariables_test_"
                                                      + std::to_string( i ) + "_"
                                                      + std::to_string( j ) + "_"
                                                      + std::to_string( k ) + ".dat", outputPath );

                // Compute difference w.r.t. benchmark using the interpolators we created
                if( generateAndCompareToBenchmark )
                {
                    std::map< double, Eigen::VectorXd> stateDifference;
                    std::map< double, Eigen::VectorXd> depVarDifference;

                    for( auto stateIterator = stateHistory.begin(); stateIterator != stateHistory.end(); stateIterator++ )
                    {
                        stateDifference[ stateIterator->first ] = stateIterator->second - benchmarkInterpolators.at( 0 )->interpolate(
                                    stateIterator->first );
                    }

                    for( auto depVarIterator = dependentVariableHistory.begin(); depVarIterator != dependentVariableHistory.end();
                         depVarIterator++ )
                    {
                        depVarDifference[ depVarIterator->first ] = depVarIterator->second - benchmarkInterpolators.at( 1 )->interpolate(
                                    depVarIterator->first );
                    }

                    // Write maps to files
                    input_output::writeDataMapToTextFile( stateDifference,
                                                          "stateDifferenceBenchmark_"
                                                          + std::to_string( i ) + "_"
                                                          + std::to_string( j ) + "_"
                                                          + std::to_string( k ) + ".dat", outputPath );

                    input_output::writeDataMapToTextFile( depVarDifference,
                                                          "dependentVariablesDifferenceBenchmark_"
                                                          + std::to_string( i ) + "_"
                                                          + std::to_string( j ) + "_"
                                                          + std::to_string( k ) + ".dat", outputPath );

                    // Write the number of function evaluations to a file for comparison of different integrators
                    unsigned int numberOfEvaluations =
                            prob.getLastRunDynamicsSimulator( )->getCumulativeNumberOfFunctionEvaluations( ).rbegin( )->second;

                    std::map< unsigned int, std::string > numberOfEvaluationsMap;
                    numberOfEvaluationsMap[numberOfEvaluations] = "";
                    input_output::writeDataMapToTextFile( numberOfEvaluationsMap,
                                                          "numberOfFunctionEvaluations.dat",
                                                          outputPath );

                }

            }
        }
    }

}
