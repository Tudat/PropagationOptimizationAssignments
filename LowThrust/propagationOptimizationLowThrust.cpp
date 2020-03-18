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

/*!
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
 */
void createSimulationSettings(
        NamedBodyMap& bodyMap,
        std::shared_ptr< MultiTypePropagatorSettings< double > >& propagatorSettings,
        std::shared_ptr< IntegratorSettings< > >& integratorSettings,
        const double specificImpulse,
        const double minimumMarsDistance,
        const double timeBuffer )
{

    // DEFINE PROBLEM INDEPENDENT VARIABLES HERE:
    std::vector< double > trajectoryParameters =
    { 1.9567061e+03,  3.8159413e+02,  0, 8.9057206e+03,  2.6738965e+03, -2.9315045e+03, 1.5350545e+03, -3.8783905e+03,  4.3249334e+03 };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            SIMULATION SETTINGS            /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define vehicle settings
    double vehicleMass = 4.0E3;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT                   //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create solar system bodies
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );

    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate );
    bodyMap = createBodies( bodySettings );

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
    ///////////////////////             RUN SIMULATION FOR VARIOUS SETTINGS            ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define translational state propagation settings
    TranslationalPropagatorType propagatorType = cowell;
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalStatePropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationSettingsMap, bodiesToPropagate, systemInitialState,
                terminationSettings, propagatorType );


    // Define full propagation settings
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList =
    { translationalStatePropagatorSettings, massPropagatorSettings };

    propagatorSettings =
            std::make_shared< MultiTypePropagatorSettings< double > >(
                propagatorSettingsList, terminationSettings, dependentVariablesToSave );
    integratorSettings =
            std::make_shared< IntegratorSettings< > >( rungeKutta4, initialPropagationTime, 3600.0 );
}

/*!
 *
 *   The trajectory of the vehicle is determined by its departure and arrival time (which define the initial and final states)
 *   as well as the free parameters of the shaping method. The free parameters of the shaping method defined here are the same
 *   as for the 'higher-order solution' in Section V.A of Gondelach and Noomen (2015). The free parameters define the amplitude
 *   of specific types of velocity shaping functions. The low-thrust hodographic trajectory is parameterized by the values of
 *   the vector trajectoryParameters
 *
 *    The entries of the vector 'trajectoryParameters' contains the following:
 *
 *   - Entry 0: Departure time (from Earth's center-of-mass) in Julian days since J2000
 *   - Entry 1: Time-of-flight from Earth's center-of-mass to Mars' center-of-mass, in Julian days
 *   - Entry 2: Number of revolutions around the Sun
 *   - Entry 3-8: Free parameters of hodographic shape (entries 3-4, 5-6 and 7-8 for radial, normal and axial directions,
 *     respectively)
 *
 */
int main( )
{
    spice_interface::loadStandardSpiceKernels( );
    std::string outputPath = tudat_applications::getOutputPath( "LowThrustDesignSpace/" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////   CREATE SIMULATION ENVIRONMENT                     ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double specificImpulse = 3000.0;
    double minimumMarsDistance = 5.0E7;
    double timeBuffer = 30.0 * physical_constants::JULIAN_DAY;

    // Create simulation environment
    NamedBodyMap bodyMap;
    std::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings;
    std::shared_ptr< IntegratorSettings< > > integratorSettings;
    createSimulationSettings(
                bodyMap, propagatorSettings, integratorSettings,
                specificImpulse, minimumMarsDistance, timeBuffer );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////   PERFORM MONTE-CARLO DESIGN SPACE EXPLORATION               //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define decision variable range
    std::vector< std::pair< double, double > > decisionVariableRange =
    {{ 0.0, 6000.0 }, { 100.0, 800.0 }, { 0, 2.999999999 }, { -10000, 10000 },
     { -10000, 10000 }, { -10000, 10000 }, { -10000, 10000 }, { -10000, 10000 }, { -10000, 10000 } };

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
    LowThrustProblem lowThrustProblem{
        bodyMap, integratorSettings, propagatorSettings, decisionVariableRange, specificImpulse, minimumMarsDistance, timeBuffer };

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
        lowThrustProblem.fitness( decisionVariables );

        std::vector< double > objectives = lowThrustProblem.getLastRunObjectives( );
        std::vector< double > constraints = lowThrustProblem.getLastRunConstraints( );

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

//    pagmo::problem prob{ lowThrustProblem };

//    // Solve using DE algorithm
//    pagmo::algorithm algo{ pagmo::de( ) };

//    // Create island with 1000 individuals
//    pagmo::island isl = pagmo::island{ algo, prob, 100 };

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

