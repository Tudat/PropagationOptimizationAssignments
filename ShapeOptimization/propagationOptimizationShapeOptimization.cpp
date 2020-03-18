/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "pagmo/island.hpp"
#include "pagmo/problem.hpp"
#include "pagmo/algorithms/de.hpp"

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


void createSimulationSettings(
        NamedBodyMap& bodyMap,
        std::shared_ptr< TranslationalStatePropagatorSettings< double > >& propagatorSettings,
        std::shared_ptr< IntegratorSettings< > >& integratorSettings,
        const double vehicleDensity )
{
    // DEFINE PROBLEM INDEPENDENT VARIABLES HERE:
    std::vector< double > shapeParameters =
    { 8.148730872315355, 2.720324489288032, 0.2270385167794302, -0.4037530896422072, 0.2781438040896319, 0.4559143679738996 };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            SIMULATION SETTINGS            /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define simulation termination settings
    double maximumDuration = 86400.0;
    double terminationAltitude = 25.0E3;

    // Vehicle properties
    // Set simulation start epoch.
    double simulationStartEpoch = 0.0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Define simulation body settings.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Moon" );

    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate );
    for( unsigned int j = 0; j < bodiesToCreate.size( ); j++ )
    {
        bodySettings[ bodiesToCreate.at( j ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
        bodySettings[ bodiesToCreate.at( j ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
    }

    // Create Earth object
    bodyMap = simulation_setup::createBodies( bodySettings );
    addCapsuleToBodyMap( bodyMap, shapeParameters, vehicleDensity );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationSettingsMap;

    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfCapsule;
    accelerationsOfCapsule[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationsOfCapsule[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
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

    // Create propagator settings
    TranslationalPropagatorType propagatorType = cowell;
    propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationSettingsMap, bodiesToPropagate, systemInitialState,
                terminationSettings, propagatorType, dependentVariablesToSave );


    // Create integrator settings
    integratorSettings =
            std::make_shared< IntegratorSettings< > >( rungeKutta4, simulationStartEpoch, 1.0 );
}

/*!
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
 */
int main()
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );
    std::string outputPath = tudat_applications::getOutputPath( "ShapeOptimizationDesignSpace/" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////   CREATE SIMULATION ENVIRONMENT                     ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double vehicleDensity = 250.0;

    NamedBodyMap bodyMap;
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings;
    std::shared_ptr< IntegratorSettings< > > integratorSettings;
    createSimulationSettings(
                bodyMap, propagatorSettings, integratorSettings, vehicleDensity );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////   PERFORM MONTE-CARLO DESIGN SPACE EXPLORATION               //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define decision variable range
    std::vector< std::pair< double, double > > decisionVariableRange =
    { { 3.5, 10.0 }, { 2.0, 3.0 }, { 0.1, 5.0 }, { -55.0 * PI / 180.0, -10.0 * PI / 180.0 }, { 0.01, 0.5 },
      { 0.0, 30.0 * PI / 180.0 } };

    // Create random number generators for decision variables
    std::vector< std::function< double( ) > > parameterMonteCarloFunctions;
    for( unsigned int i = 0; i < decisionVariableRange.size( ); i++ )
    {
        parameterMonteCarloFunctions.push_back(
                    statistics::createBoostContinuousRandomVariableGeneratorFunction(
                        statistics::uniform_boost_distribution,
        { decisionVariableRange.at( i ).first, decisionVariableRange.at( i ).second }, i ) );
    }

    // Create problem class
    ShapeOptimizationProblem shapeOptimizationProblem{
        bodyMap, integratorSettings, propagatorSettings, decisionVariableRange, vehicleDensity };

    // Run Monte Carlo design space exploration
    for( int i = 0; i < 10000; i++ )
    {
        std::cout<<i<<std::endl;
        std::vector< double > decisionVariables;
        for( unsigned int j = 0; j < parameterMonteCarloFunctions.size( ); j++ )
        {
            decisionVariables.push_back( parameterMonteCarloFunctions.at( j )( ) );
        }

        // Construct problem and propagate trajectory using defined settings
        shapeOptimizationProblem.fitness( decisionVariables );

        std::vector< double > objectives = shapeOptimizationProblem.getLastRunObjectives( );
        std::vector< double > constraints = shapeOptimizationProblem.getLastRunConstraints( );

        input_output::writeMatrixToFile(
                    utilities::convertStlVectorToEigenVector( decisionVariables ),
                    "decisionVariables_" + std::to_string( i ) + ".dat", 16, outputPath );
        input_output::writeMatrixToFile(
                    utilities::convertStlVectorToEigenVector( objectives ),
                    "objectives_" + std::to_string( i ) + ".dat", 16, outputPath );
        input_output::writeMatrixToFile(
                    utilities::convertStlVectorToEigenVector( constraints ),
                    "constraints_" + std::to_string( i ) + ".dat", 16, outputPath );
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////   EXAMPLE CODE FOR A SIMPLE OPTIMIZATION            ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//    pagmo::problem prob{ shapeOptimizationProblem };

//    // Solve using DE algorithm
//    pagmo::algorithm algo{ pagmo::de( ) };

//    // Create island with 1000 individuals
//    pagmo::island isl = pagmo::island{ algo, prob, 1000 };

//    // Evolve for 1000 generations
//    for( int i = 1; i <= 100; i++ )
//    {
//        isl.evolve( );
//        while( isl.status( ) != pagmo::evolve_status::idle &&
//               isl.status( ) != pagmo::evolve_status::idle_error )
//        {
//            isl.wait( );
//        }
//        isl.wait_check( ); // Raises errors
//        // Print current optimum to console
//        std::cout << "Minimum: " <<i<<" "<<std::setprecision( 16 ) <<"f= "<< isl.get_population().champion_f()[0] <<", x="<<
//                     isl.get_population().champion_x()[0] <<" y="<<isl.get_population().champion_x()[1] <<std::endl;

//    }

}

