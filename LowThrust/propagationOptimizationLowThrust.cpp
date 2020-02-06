/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "lowThrust.h"

using namespace tudat_applications::PropagationOptimization2020;

//! Function to write properties of shape-based hodographic solution to file
/*!
 *
 *  Function to write properties of shape-based hodographic solution to file. For a given trajectory shape, this function writes:
 *
 *   - hodographicTrajectory.dat: Cartesian states of semi-analytical trajectory
 *   - hodographicThrustAcceleration.dat: Thrust acceleration in inertial, Cartesian, coordinates, along the semi-analytical
 *                                        trajectory.
 *
 *  NOTE: The independent variable (first column) does not represent teh usual time (seconds since J2000), but instead denotes
 *  the time since departure.
 *
 *  These files are written to the directory specified by the 'outputPath' variables
 */
void printSemiAnalyticalHodographicShapeToFile(
        const std::shared_ptr< HodographicShaping > shapeObject,
        const std::vector< double >& trajectoryParameters,
        const double specificImpulse,
        const std::string outputPath )
{
    std::vector< double > epochs;
    double startTime = 0.0;
    double finalTime = getTrajectoryTimeOfFlight( trajectoryParameters );
    int numberOfDataPoints = 10000;
    double stepSize = ( finalTime - startTime ) / static_cast< double >( numberOfDataPoints - 1 );
    for( int i = 0; i < numberOfDataPoints; i++ )
    {
        epochs.push_back( startTime + static_cast< double >( i ) * stepSize );
    }

    std::map< double, Eigen::VectorXd > thrustAccelerationProfile;
    shapeObject->getThrustAccelerationProfile(
                epochs, thrustAccelerationProfile, [=](const double){return specificImpulse; }, nullptr );

    std::map< double, Eigen::Vector6d > propagatedTrajectory;
    shapeObject->getTrajectory( epochs, propagatedTrajectory );

    input_output::writeDataMapToTextFile(
                propagatedTrajectory, "hodographicTrajectory.dat", outputPath );
    input_output::writeDataMapToTextFile(
                thrustAccelerationProfile, "hodographicThrustAcceleration.dat", outputPath );
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
 *      if j=4      : a fixed-step-size RK4 integrator is used, with step-size 2 hours*2^(k)
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
        double timeStep = 7200.0 * std::pow( 2, k );
        return std::make_shared< IntegratorSettings< > >( rungeKutta4, simulationStartEpoch, timeStep );
    }
}

//! Function to retrieve the dependent variable save settings for the current simulation.
/*!
 * This function returns a shared pointer to a DependentVariableSaveSettings object, containing the save settings
 * to save the altitude, relative speed (w.r.t. Moon center of mass) and flight path angle of the vehicle
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
                                          relative_distance_dependent_variable, "Vehicle", "Earth" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          relative_distance_dependent_variable, "Vehicle", "Sun" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          relative_distance_dependent_variable, "Vehicle", "Mars" ) );
    return std::make_shared< DependentVariableSaveSettings >( dependentVariablesList, false );
}

//! Function to generate to accurate benchmarks.
/*!
 * This function runs two propagations with two different integrator settings that serve as benchmarks for
 * the nominal runs. To be able to compare these, the function returns the two interpolators pertaining
 * to the state and dependent variables of one of the benchmarks. The states are written to files, as well
 * as the difference in state and dependent variables between the two benchmarks.
 *
 *  The following files are written to files by this function (to the directory LowThrust/benchmarks/...:
 *
 *  - benchmarkStates_1.dat, benchmarkStates_2.dat The numerically propagated states from the two propagations
 *  - benchmarkDependent_1.dat, benchmarkDependent_2.dat The dependent variables from the two propagations
 *  - benchmarkStateDifference.dat Difference between the Cartesian states of the two benchmark runs
 *  - benchmarkDependentDifference.dat  Difference between the dependent variables of the two benchmark runs
 *
 *
 *  CODING NOTE: THIS FUNCTION CAN BE EXTENDED TO GENERATE A MORE ROBUST BENCHMARK (USING MORE THAN 2 RUNS)
 *
 * \param simulationStartEpoch The start time of the simulation in seconds.
 * \param bodyMap NamedBodyMap containing the bodies in the simulation.
 * \param benchmarkPropagatorSettings Shared pointer to a translational propagator settings object,
 * which is used to run the benchmark propagations.
 * \param trajectoryParameters The vector of doubles that represents the trajectory parameters for the capsule.
 * \param outputPath String containing the path to the output directory.
 * \return
 */
std::vector< std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > > generateBenchmarks(
        const double simulationStartEpoch, const double specificImpulse, const double minimumMarsDistance,
        const double timeBuffer,const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< MultiTypePropagatorSettings< double > > benchmarkPropagatorSettings,
        std::vector< double > trajectoryParameters, std::string outputPath )
{
    // Create integrator settings for 1st run
    double firstBenchmarkStepSize = 86400.0;
    std::shared_ptr< IntegratorSettings< > > benchmarkIntegratorSettings;
    benchmarkIntegratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                simulationStartEpoch, firstBenchmarkStepSize, RungeKuttaCoefficients::rungeKutta87DormandPrince,
                firstBenchmarkStepSize, firstBenchmarkStepSize,
                std::numeric_limits< double >::infinity( ), std::numeric_limits< double >::infinity( ) );

    LowThrustProblem probBenchmarkFirst{
        bodyMap, benchmarkIntegratorSettings, benchmarkPropagatorSettings, specificImpulse, minimumMarsDistance, timeBuffer };

    std::cout << "Running first benchmark..." << std::endl;
    probBenchmarkFirst.fitness( trajectoryParameters );

    // Create integrator settings for 2nd run
    double secondBenchmarkStepSize = 2.0 * 86400.0;
    benchmarkIntegratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                simulationStartEpoch, secondBenchmarkStepSize, RungeKuttaCoefficients::rungeKutta87DormandPrince,
                secondBenchmarkStepSize, secondBenchmarkStepSize,
                std::numeric_limits< double >::infinity( ), std::numeric_limits< double >::infinity( ) );
    LowThrustProblem probBenchmarkSecond{
        bodyMap, benchmarkIntegratorSettings, benchmarkPropagatorSettings, specificImpulse, minimumMarsDistance, timeBuffer  };

    std::cout << "Running second benchmark..." << std::endl;
    probBenchmarkSecond.fitness( trajectoryParameters );

    // Retrieve states and dependent variables for both runs
    std::map< double, Eigen::VectorXd > firstBenchmarkStates = probBenchmarkFirst.getLastRunPropagatedStateHistory( );
    std::map< double, Eigen::VectorXd > secondBenchmarkStates = probBenchmarkSecond.getLastRunPropagatedStateHistory( );


    std::map< double, Eigen::VectorXd > firstBenchmarkDependent = probBenchmarkFirst.getLastRunDependentVariableHistory( );
    std::map< double, Eigen::VectorXd > secondBenchmarkDependent = probBenchmarkSecond.getLastRunDependentVariableHistory( );
    \
    // Put the benchmark data in a separate directory
    outputPath.append("/benchmarks/");

    // Write the state maps of both benchmarks to files
    input_output::writeDataMapToTextFile( firstBenchmarkStates, "benchmark1.dat", outputPath );
    input_output::writeDataMapToTextFile( secondBenchmarkStates, "benchmark2.dat", outputPath );

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
 *   This function computes the dynamics of an interplanetary low-thrust trajectory, using a thrust profile determined from
 *   a Hodographic shaping method (see Gondelach and Noomen, 2015). This file propagates the dynamics using a variety of integrator
 *   and propagator settings (see comments under "RUN SIMULATION FOR VARIOUS SETTINGS"). For each run, the differences w.r.t. a
 *   benchmark propagation are computed, providing a proxy for setting quality.
 *
 *   The low-thrust trajectory computed by the shape-based method starts at the Earth's center of mass, and terminates at Mars'
 *   center of mass.
 *
 *   The vehicle starts on the Hodographic low-thrust trajectory, 30 days (defined by the timeBuffer variable) after it 'departs'
 *   the Earth's center of mass.
 *
 *   The propagation is terminated as soon as one of the following conditions is met (see
 *   getPropagationTerminationSettings function):
 *
 *   - Distance to Mars < 50000 km
 *   - Propagation time > Time-of-flight of hodographic trajectory
 *
 *   This propagation assumes only point mass gravity by the Sun and thrust acceleration of the vehicle
 *   (see block 'CREATE ACCELERATIONS'). Both the translational dynamics and mass of the vehicle are propagated,
 *   using a fixed specific impulse.
 *
 *   The trajectory of the capsule is determined by its departure and arrival time (which define the initial and final states)
 *   as well as the free parameters of the shaping method. The free parameters of the shaping method defined here are the same
 *   as for the 'higher-order solution' in Section V.A of Gondelach and Noomen (2015). The free parameters define the amplitude
 *   of specific types of velocity shaping functions. The low-thrust hodographic trajectory is parameterized by the values of
 *   the vector trajectoryParameters
 *
 *    The entries of the vector 'trajectoryParameters' contains the following:
 *
 *   - Entry 0: Departure time (from Earth's center-of-mass) in Julian days since J2000
 *   - Entry 1: Time-of-flight from Earth's center-of-mass to Mars' center-of-mass, in Julian days
 *   - Entry 2-6: Thrust angle theta, at nodes 1-5 (in order)
 *
 *   Details on the outputs written by this file can be found:
 *
 *      Benchmark data: comments for 'generateBenchmarks' function
 *      Results for integrator/propagator variations: comments under "RUN SIMULATION FOR VARIOUS SETTINGS"
 *      Trajectory for semi-analytical hodographic shape-based solution: Comments with, and call to,
 *                                                                       printSemiAnalyticalHodographicShapeToFile function
 *
 */
int main( )
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // DEFINE PROBLEM INDEPENDENT VARIABLES HERE:
    std::vector< double > trajectoryParameters =
    { 1.9567061e+03,  3.8159413e+02,  0, 8.9057206e+03,  2.6738965e+03, -2.9315045e+03, 1.5350545e+03, -3.8783905e+03,  4.3249334e+03 };

    std::string outputPath = tudat_applications::getOutputPath( "LowThrust" );
    bool generateAndCompareToBenchmark = true;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            SIMULATION SETTINGS            /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define vehicle settings
    double vehicleMass = 4.0E3;
    double specificImpulse = 3000.0;
    double minimumMarsDistance = 5.0E7;
    double timeBuffer = 30.0 * physical_constants::JULIAN_DAY;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT                   //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create solar system bodies
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Sun" );

    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate );
    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationSettingsMap;

    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );

    accelerationsOfVehicle[ "Vehicle" ].push_back(
                getThrustAccelerationSettingsFromParameters( trajectoryParameters, bodyMap ) );
    accelerationSettingsMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Define propagator settings variables.
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Sun" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////   RETRIEVE DATA FOR PROPAGATION SETTINGS            ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double initialPropagationTime = getTrajectoryInitialTime( trajectoryParameters, timeBuffer );
    std::shared_ptr< PropagationTerminationSettings > terminationSettings = getPropagationTerminationSettings(
                trajectoryParameters, minimumMarsDistance, 0.0 );
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave = getDependentVariableSaveSettings();
    Eigen::Vector6d systemInitialState = getHodographicLowThrustStateAtEpoch(
                trajectoryParameters, bodyMap, initialPropagationTime );

    // Create propagator settings for mass (constant for all simulations)
    simulation_setup::SelectedMassRateModelMap massRateModelSettings;
    massRateModelSettings[ "Vehicle" ].push_back( std::make_shared< FromThrustMassModelSettings >( 1 ) );
    std::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings =
            std::make_shared< MassPropagatorSettings< double > >(
                std::vector< std::string >{ "Vehicle" }, massRateModelSettings,
                ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ), terminationSettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////  IF DESIRED, GENERATE BENCHMARK                            ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector< std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > > benchmarkInterpolators;
    if( generateAndCompareToBenchmark )
    {
        // Define translational state propagation settings
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalStatePropagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies, accelerationSettingsMap, bodiesToPropagate, systemInitialState,
                    terminationSettings, cowell );

        // Define full propagation settings
        std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList =
        { translationalStatePropagatorSettings, massPropagatorSettings };
        std::shared_ptr< MultiTypePropagatorSettings< double > > benchmarkPropagatorSettings =
                std::make_shared< MultiTypePropagatorSettings< double > >(
                    propagatorSettingsList, terminationSettings, dependentVariablesToSave );

        // Generate benchmark data
        benchmarkInterpolators = generateBenchmarks(
                    initialPropagationTime, specificImpulse, minimumMarsDistance, timeBuffer, bodyMap, benchmarkPropagatorSettings,
                    trajectoryParameters, outputPath );

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     PRINT DATA FOR HODOGRAPHIC SEMI-ANALYTICAL METHOD         /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    LowThrustProblem semiAnalyticalProblem{
        bodyMap, nullptr, nullptr, specificImpulse, minimumMarsDistance, timeBuffer, false };
    semiAnalyticalProblem.fitness( trajectoryParameters );
    std::shared_ptr< HodographicShaping > shapeObject = semiAnalyticalProblem.getHodographicShaping( );

    printSemiAnalyticalHodographicShapeToFile(
                shapeObject, trajectoryParameters, specificImpulse, outputPath + "HodographicSemiAnalytical" );

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
    //!     SimulationOutput/LowThrust/prop_i/int_j/setting_k/
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
        // Define translational state propagation settings
        TranslationalPropagatorType propagatorType = propagatorTypes.at( i );
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalStatePropagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies, accelerationSettingsMap, bodiesToPropagate, systemInitialState,
                    terminationSettings, propagatorType );

        // Define full propagation settings
        std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList =
        { translationalStatePropagatorSettings, massPropagatorSettings };
        std::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings =
                std::make_shared< MultiTypePropagatorSettings< double > >(
                    propagatorSettingsList, terminationSettings, dependentVariablesToSave );

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
                            "LowThrust/prop_" + std::to_string( i ) + "/int_" + std::to_string( j ) + "/setting_" + std::to_string( k ) + "/" );

                // Create integrator settings
                std::shared_ptr< IntegratorSettings< > > integratorSettings = getIntegratorSettings( i, j, k, initialPropagationTime );

                // Construct problem and propagate trajectory using defined settings
                LowThrustProblem prob{
                    bodyMap, integratorSettings, propagatorSettings, specificImpulse, minimumMarsDistance, timeBuffer };
                prob.fitness( trajectoryParameters );

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
