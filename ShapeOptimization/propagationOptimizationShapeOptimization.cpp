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

//! Get the propagation termination settings for the state propagation
/*!
 * This function returns a shared pointer to a PropagationTerminationSettings object, containing settings termination on:
 *
 *      altitude                (<terminationAltitude)
 *      total propagation time  (>maximumDuration)
 *
 * The settings are such that the propagation terminates once at least one of these conditions has been met
 * \param initialTime Start time of the simulation in seconds.
 * \param maximumDuration Time in seconds, specifying the maximum time duration before which the
 * simulation should stop.
 * \param terminationAltitude Altitude in meters, specifying the maximum altitude before which the
 * simulation should stop.
 * \return Shared pointer to the PropagationTerminationSettings object.
 */
std::shared_ptr< PropagationTerminationSettings > getPropagationTerminationSettings(
        const double initialTime,
        const double maximumDuration,
        const double terminationAltitude )
{
    // Define termination conditions
    std::vector< std::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;

    // Add altitude termination condition
    std::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
            std::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "Capsule", "Earth" );
    terminationSettingsList.push_back(
                std::make_shared< PropagationDependentVariableTerminationSettings >( terminationDependentVariable, terminationAltitude, true ) );

    // Add time termination condision
    terminationSettingsList.push_back( std::make_shared< PropagationTimeTerminationSettings >( initialTime + maximumDuration ) );

    return std::make_shared< PropagationHybridTerminationSettings >( terminationSettingsList, true );
}

//! Function to retrieve the integrator settings for the current run.
/*!
 * This function returns a shared pointer to an IntegratorSettings object, based on the indices passed to
 * the function. The first index is used to determine which type is selected from the vector of integrator
 * types which is passed as the last parameter. The seconds index determines what tolerance is used for the
 * variable step size integrators.
 *
 * The code, as provided, runs the following:
 *      if j=0,1,2,3: a variable-step-size, multi-stage integrator is used (see multiStageTypes list for specific type),
 *                    with tolerances 10^(-10+k)
 *      if j=4      : a fixed-step-size RK4 integrator is used, with step-size 10 * 2^(k)
 *
 * CODING NOTE: THIS FUNCTION SHOULD BE EXTENDED TO USE MORE INTEGRATORS FOR ASSIGNMENT 1
 *
 * \param i Index specifying which kind of propagator is used
 * \param j Index specifying which kind of integrator is used (see above)
 * \param k Index that is used to specify different tolerances for the same integrator (see above)
 * \param simulationStartEpoch The start time of the simulation in seconds.
 * \return Shared pointer to the IntegratorSettings object.
 */
std::shared_ptr< IntegratorSettings< > > getIntegratorSettings(
        unsigned int i, unsigned int j, unsigned int k, double simulationStartEpoch )
{
    // Define list of multi-stage integrators (for convenience)
    std::vector< RungeKuttaCoefficients::CoefficientSets > multiStageTypes =
    { RungeKuttaCoefficients::rungeKuttaFehlberg45,
      RungeKuttaCoefficients::rungeKuttaFehlberg56,
      RungeKuttaCoefficients::rungeKuttaFehlberg78,
      RungeKuttaCoefficients::rungeKutta87DormandPrince };

    // Create integrator settings (multi-stage variable-step)
    if( j < 4 )
    {
        // Extract integrator type and tolerance for current run
        RungeKuttaCoefficients::CoefficientSets currentCoefficientSet = multiStageTypes.at( j );
        double currentTolerance = std::pow( 10.0, ( -10.0 + static_cast< double >( k ) ) );

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
        double timeStep = 10.0 * std::pow( 2, k );
        return std::make_shared< IntegratorSettings< > >( rungeKutta4, simulationStartEpoch, timeStep );
    }
}

//! Function to retrieve the dependent variable save settings for the current simulation.
/*!
 * This function returns a shared pointer to a DependentVariableSaveSettings object, containing the save settings
 * to save the Mach number and altitude of the capsule.
 *
 *  CODING NOTE: THIS FUNCTION SHOULD BE EXTENDED TO SAVE MORE DEPENDENT VARIABLES
 *
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
 *  This function runs two propagations with two different integrator settings that serve as benchmarks for
 *  the nominal runs. To be able to compare these, the function returns the two interpolators pertaining
 *  to the state and dependent variables of one of the benchmarks. The states are written to files, as well
 *  as the difference in state and dependent variables between the two benchmarks.
 *
 *  The following files are written to files by this function (to the directory ShapeOptimziation/benchmarks/...:
 *
 *  - benchmarkStates_1.dat, benchmarkStates_2.dat The numerically propagated states from the two propagations
 *  - benchmarkDependent_1.dat, benchmarkDependent_2.dat The dependent variables from the two propagations
 *  - benchmarkStateDifference.dat Difference between the Cartesian states of the two benchmark runs
 *  - benchmarkDependentDifference.dat  Difference between the dependent variables of the two benchmark runs
 *
 *
 *
 *  CODING NOTE: THIS FUNCTION CAN BE EXTENDED TO GENERATE A MORE ROBUST BENCHMARK (USING MORE THAN 2 RUNS)
 *
 * \param simulationStartEpoch The start time of the simulation in seconds.
 * \param bodyMap NamedBodyMap containing the bodies in the simulation.
 * \param benchmarkPropagatorSettings Shared pointer to a translational propagator settings object,
 * which is used to run the benchmark propagations.
 * \param shapeParameters The vector of doubles that represents the shape parameters for the capsule.
 * \param outputPath String containing the path to the output directory.
 * \return Interpolators providing values for state and dependent variables of the benchmark run
 */
std::vector< std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > > generateBenchmarks(
        double simulationStartEpoch, const double vehicleDensity, simulation_setup::NamedBodyMap bodyMap,
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > benchmarkPropagatorSettings,
        std::vector< double > shapeParameters, std::string outputPath )
{
    // Create integrator settings for 1st run
    double firstBenchmarkStepSize = 2.0;
    std::shared_ptr< IntegratorSettings< > > benchmarkIntegratorSettings;
    benchmarkIntegratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                simulationStartEpoch, firstBenchmarkStepSize, RungeKuttaCoefficients::rungeKutta87DormandPrince,
                firstBenchmarkStepSize, firstBenchmarkStepSize,
                std::numeric_limits< double >::infinity( ), std::numeric_limits< double >::infinity( ) );

    ShapeOptimizationProblem probBenchmarkFirst{ bodyMap, benchmarkIntegratorSettings, benchmarkPropagatorSettings, vehicleDensity };

    std::cout << "Running first benchmark..." << std::endl;
    probBenchmarkFirst.fitness( shapeParameters );

    // Create integrator settings for 2nd run
    double secondBenchmarkStepSize = 4.0;
    benchmarkIntegratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                simulationStartEpoch, secondBenchmarkStepSize, RungeKuttaCoefficients::rungeKutta87DormandPrince,
                secondBenchmarkStepSize, secondBenchmarkStepSize,
                std::numeric_limits< double >::infinity( ), std::numeric_limits< double >::infinity( ) );
    ShapeOptimizationProblem probBenchmarkSecond{ bodyMap, benchmarkIntegratorSettings, benchmarkPropagatorSettings, vehicleDensity };

    std::cout << "Running second benchmark..." << std::endl;
    probBenchmarkSecond.fitness( shapeParameters );

    // Retrieve states and dependent variables for both runs
    std::map< double, Eigen::VectorXd > firstBenchmarkStates = probBenchmarkFirst.getLastRunPropagatedStateHistory( );
    std::map< double, Eigen::VectorXd > secondBenchmarkStates = probBenchmarkSecond.getLastRunPropagatedStateHistory( );

    std::map< double, Eigen::VectorXd > firstBenchmarkDependent = probBenchmarkFirst.getLastRunDependentVariableHistory( );
    std::map< double, Eigen::VectorXd > secondBenchmarkDependent = probBenchmarkSecond.getLastRunDependentVariableHistory( );
    \
    // Put the benchmark data in a separate directory
    outputPath.append("/benchmarks/");

    // Write the state and dependent variable maps of both benchmarks to files
    input_output::writeDataMapToTextFile( firstBenchmarkStates, "benchmarkStates_1.dat", outputPath );
    input_output::writeDataMapToTextFile( secondBenchmarkStates, "benchmarkStates_2.dat", outputPath );
    input_output::writeDataMapToTextFile( firstBenchmarkDependent, "benchmarkDependent_1.dat", outputPath );
    input_output::writeDataMapToTextFile( secondBenchmarkDependent, "benchmarkDependent_2.dat", outputPath );

    // Create 8th-order Lagrange interpolators for states and dependent variables of both runs
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
    for( auto iterator = secondBenchmarkStates.begin(); iterator != secondBenchmarkStates.end(); iterator++ )
    {
        statesDifference[ iterator->first ] = firstStatesInterpolator->interpolate( iterator->first ) - iterator->second;
    }

    for( auto iterator = secondBenchmarkDependent.begin(); iterator != secondBenchmarkDependent.end(); iterator++ )
    {
        dependentVariablesDifference[ iterator->first ] = firstDependentInterpolator->interpolate( iterator->first ) - iterator->second;
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

/*!
 *   This function computes the dynamics of a capsule re-entering the atmosphere of the Earth, using a variety of integrator
 *   and propagator settings (see comments under "RUN SIMULATION FOR VARIOUS SETTINGS"). For each run, the differences w.r.t. a
 *   benchmark propagation are computed, providing a proxy for setting quality.
 *
 *   The vehicle starts 120 km above the surface of the planet, with a speed of 7.83 km/s in an Earth-fixed frame (see
 *   getInitialState function).
 *
 *   The propagation is terminated as soon as one of the following conditions is met (see
 *   getPropagationTerminationSettings function):
 *
 *   - Altitude < 25 km
 *   - Propagation time > 24 hr
 *
 *   This propagation assumes only point mass gravity by the Earth and aerodynamic accelerations (see block 'CREATE ACCELERATIONS')
 *
 *   The trajectory of the capsule is heavily dependent on the shape and orientation of the vehicle. Here, the shape is
 *   determined here by the five parameters, which are used to compute the aerodynamic accelerations on the
 *   vehicle using a modified Newtonian flow (see also Dirkx and Mooij, 2018). The bank angle and sideslip angles are set to zero.
 *   The vehicle shape and angle of attack are defined by values in the vector shapeParameters.
 *
 *   The entries of the vector 'shapeParameters' contains the following:
 *
 *   - Entry 0:  Nose radius
 *   - Entry 1:  Middle radius
 *   - Entry 2:  Rear length
 *   - Entry 3:  Rear angle
 *   - Entry 4:  Side radius
 *   - Entry 5:  Constant Angle of Attack
 *
 *   Details on the outputs written by this file can be found:
 *
 *      Benchmark data: comments for 'generateBenchmarks' function
 *      Results for integrator/propagator variations: comments under "RUN SIMULATION FOR VARIOUS SETTINGS"
 *
 */
int main()
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // DEFINE PROBLEM INDEPENDENT VARIABLES HERE:
    std::vector< double > shapeParameters =
    { 8.148730872315355, 2.720324489288032, 0.2270385167794302, -0.4037530896422072, 0.2781438040896319, 0.4559143679738996 };

    std::string outputPath = tudat_applications::getOutputPath( "ShapeOptimization" );
    bool generateAndCompareToBenchmark = true;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            SIMULATION SETTINGS            /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define simulation termination settings
    double maximumDuration = 86400.0;
    double terminationAltitude = 25.0E3;

    // Vehicle properties
    double vehicleDensity = 250.0;

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

    //! Create capsule, and add to body map
    //! NOTE: When making any modifications to the capsule vehicle, do NOT make them in this main function,
    //! but in the addCapsuleToBodyMap function
    addCapsuleToBodyMap( bodyMap, shapeParameters, vehicleDensity );

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

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////   RETRIEVE DATA FOR PROPAGATION SETTINGS            ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::shared_ptr< PropagationTerminationSettings > terminationSettings = getPropagationTerminationSettings(
                simulationStartEpoch, maximumDuration, terminationAltitude );
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave = getDependentVariableSaveSettings();
    Eigen::Vector6d systemInitialState = getInitialState( simulationStartEpoch, bodyMap );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////  IF DESIRED, GENERATE BENCHMARK                            ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector< std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > > benchmarkInterpolators;
    if( generateAndCompareToBenchmark )
    {
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > benchmarkPropagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies, accelerationSettingsMap, bodiesToPropagate, systemInitialState,
                    terminationSettings, cowell, dependentVariablesToSave );
        benchmarkInterpolators = generateBenchmarks(simulationStartEpoch, vehicleDensity, bodyMap, benchmarkPropagatorSettings,
                                                    shapeParameters, outputPath );
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             RUN SIMULATION FOR VARIOUS SETTINGS            ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define list of propagators (for convenience)
    std::vector< TranslationalPropagatorType > propagatorTypes =
    { cowell, encke, gauss_keplerian, gauss_modified_equinoctial,
      unified_state_model_quaternions, unified_state_model_modified_rodrigues_parameters,
      unified_state_model_exponential_map };

    //!  Code below propagates states using each propagator (index i=0..6), four multi-stage variable step-size integrators
    //!  (index j=0..3) and an RK4 integrator (j=4). For the variable-step integrators, 4 different tolerances are used (k=0..3).
    //!  For the RK4, 6 different step sizes are used (k=0..5), see use of numberOfIntegratorStepSizeSettings variable. See
    //!  getIntegratorSettings function for more details.
    //!
    //!  For each combination of i,j and k, results are written to directory:
    //!
    //!     SimulationOutput/ShapeOptimization/prop_i/int_j/setting_k/
    //!
    //!  Specifically:
    //!
    //!      stateHistory.dat                           Cartesian states as function of time
    //!      dependentVariables.dat                     Dependent variables as function of time
    //!      stateDifferenceBenchmark.dat               Difference of dependent variables w.r.t. benchmark
    //!      dependentVariablesDifferenceBenchmark.dat  Difference of states w.r.t. benchmark
    //!      numberOfFunctionEvaluations.dat            Number of function evaluations performed by propagation
    //!      propagationSuccesfull.dat                  Boolean denoting whether the propagation was succesful (if false, disregard propagation)
    //!
    //!
    //!  CODING NOTE: THE NUMBER, TYPES, SETTINGS OF PROPAGATORS/INTEGRATORS/INTEGRATOR STEPS,TOLERANCES,ETC. SHOULD BE MODIFIED FOR ASSIGNMENT 1
    //!
    unsigned int numberOfPropagators = 7;
    unsigned int numberOfIntegrators = 5;
    unsigned int numberOfIntegratorStepSizeSettings = 4;
    for( unsigned int i = 0; i < numberOfPropagators; i++ )
    {
        // Create propagator settings
        TranslationalPropagatorType propagatorType = propagatorTypes.at( i );
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies, accelerationSettingsMap, bodiesToPropagate, systemInitialState,
                    terminationSettings, propagatorType, dependentVariablesToSave );

        // Iterate over all types of integrators
        for( unsigned int j = 0; j < numberOfIntegrators; j++ )
        {
            // Change number of step sizes used for RK4
            if( j >= 4 )
            { numberOfIntegratorStepSizeSettings = 6; }
            else
            { numberOfIntegratorStepSizeSettings = 4; }

            // Iterate over all tolerances/step sizes
            for( unsigned int k = 0; k < numberOfIntegratorStepSizeSettings; k++ )
            {
                // Print status
                std::cout<<"Current run "<<i<<" "<<j<<" "<<k<<std::endl;
                outputPath = tudat_applications::getOutputPath(
                            "ShapeOptimization/prop_" + std::to_string( i ) + "/int_" + std::to_string( j ) + "/setting_" + std::to_string( k ) + "/" );

                // Create integrator settings
                std::shared_ptr< IntegratorSettings< > > integratorSettings = getIntegratorSettings( i, j, k, simulationStartEpoch );

                // Construct problem and propagate trajectory using defined settings
                ShapeOptimizationProblem prob{ bodyMap, integratorSettings, propagatorSettings, vehicleDensity };
                prob.fitness( shapeParameters );

                // Save state and dependent variable results to file
                std::map< double, Eigen::VectorXd> stateHistory = prob.getLastRunPropagatedStateHistory( );
                std::map< double, Eigen::VectorXd> dependentVariableHistory = prob.getLastRunDependentVariableHistory( );
                input_output::writeDataMapToTextFile( stateHistory,  "stateHistory.dat", outputPath );
                input_output::writeDataMapToTextFile( dependentVariableHistory, "dependentVariables.dat", outputPath );

                // Write the number of function evaluations to a file for comparison of different integrators
                int numberOfEvaluations =
                        prob.getLastRunDynamicsSimulator( )->getCumulativeNumberOfFunctionEvaluations( ).rbegin( )->second;
                input_output::writeMatrixToFile( ( Eigen::MatrixXd( 1, 1 ) << numberOfEvaluations ).finished( ),
                                                 "numberOfFunctionEvaluations.dat", 16, outputPath );

                // Write to file whether the simulation was run succesfully (true or false)
                bool succesfullyRun = prob.getLastRunDynamicsSimulator( )->integrationCompletedSuccessfully( );
                input_output::writeMatrixToFile( ( Eigen::MatrixXd( 1, 1 ) << succesfullyRun ).finished( ),
                                                 "propagationSuccesfull.dat", 16, outputPath );

                // Compare to benchmark, and write differences to files
                if( generateAndCompareToBenchmark )
                {
                    std::map< double, Eigen::VectorXd> stateDifference;
                    std::map< double, Eigen::VectorXd> depVarDifference;

                    // Compute difference w.r.t. benchmark using the interpolators we created
                    for( auto stateIterator = stateHistory.begin(); stateIterator != stateHistory.end(); stateIterator++ )
                    {
                        if( dependentVariableHistory.count( stateIterator->first ) != 0 )
                        {
                            stateDifference[ stateIterator->first ] =
                                    stateIterator->second -
                                    benchmarkInterpolators.at( 0 )->interpolate( stateIterator->first );
                            depVarDifference[ stateIterator->first ] =
                                    dependentVariableHistory.at( stateIterator->first ) -
                                    benchmarkInterpolators.at( 1 )->interpolate( stateIterator->first );
                        }
                    }

                    // Write differences w.r.t. benchmarks to files
                    input_output::writeDataMapToTextFile(
                                stateDifference, "stateDifferenceBenchmark.dat", outputPath );
                    input_output::writeDataMapToTextFile(
                                depVarDifference, "dependentVariablesDifferenceBenchmark.dat", outputPath );
                }
            }
        }
    }
}
