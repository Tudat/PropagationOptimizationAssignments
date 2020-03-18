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

#include "lunarAscent.h"

using namespace tudat_applications::PropagationOptimization2020;

//! Function to retrieve the initial Cartesian state of the vehicle.
/*!
 * Function to retrieve the initial Cartesian state of the vehicle. The spherical orbital parameters are
 * first converted to Cartesian coordinates and subsequently transformed to the global frame of reference.
 *
 * \param simulationStartEpoch The start time of the simulation in seconds.
 * \param bodyMap NamedBodyMap containing the bodies in the simulation.
 * \return Eigen Vector6d containing the system's initial state in Cartesian coordinates.
 */
Eigen::Vector6d getInitialState( double simulationStartEpoch, simulation_setup::NamedBodyMap bodyMap )
{
    // Define initial spherical elements for vehicle.
    Eigen::Vector6d ascentVehicleSphericalEntryState;
    ascentVehicleSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
            spice_interface::getAverageRadius( "Moon" ) + 100.0;
    ascentVehicleSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) =
            unit_conversions::convertDegreesToRadians( 0.6875 );
    ascentVehicleSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) =
            unit_conversions::convertDegreesToRadians( 23.4333 );
    ascentVehicleSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 10.0;
    ascentVehicleSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
            unit_conversions::convertDegreesToRadians( 90.0 );
    ascentVehicleSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) =
            unit_conversions::convertDegreesToRadians( 90.0 );

    // Convert vehicle state from spherical elements to body-fixed Cartesian elements.
    Eigen::Vector6d bodyFixedSystemInitialState = convertSphericalOrbitalToCartesianState(
                ascentVehicleSphericalEntryState );

    // Convert the state to the global (inertial) frame.
    std::shared_ptr< ephemerides::RotationalEphemeris > moonRotationalEphemeris =
            bodyMap.at( "Moon" )->getRotationalEphemeris( );
    return transformStateToGlobalFrame(
                bodyFixedSystemInitialState, simulationStartEpoch, moonRotationalEphemeris );

}

//! Get the propagation termination settings for the state propagation
/*!
 * This function returns a shared pointer to a PropagationTerminationSettings object, containing settings termination on:
 *
 *      altitude                (>terminationAltitude)
 *      altitude                (<0)
 *      total propagation time  (>maximumDuration)
 *      vehicle mass            (<vehicleDryMass)
 *
 * The settings are such that the propagation terminates once at least one of these conditions has been met
 * \param simulationStartEpoch Start time of the simulation in seconds.
 * \param maximumDuration Time in seconds, specifying the maximum time duration before which the
 * simulation should stop.
 * \param terminationAltitude Altitude in meters, specifying the maximum altitude before which the
 * simulation should stop.
 * \param vehicleDryMass Dry mass of the vehicle in kg. This is value is used to create a termination
 * condition that mandates the simulation to stop once all fuel has been used up.
 * \return Shared pointer to the PropagationTerminationSettings object.
 */
std::shared_ptr< PropagationTerminationSettings > getPropagationTerminationSettings( const double simulationStartEpoch,
                                                                                     const double maximumDuration,
                                                                                     const double terminationAltitude,
                                                                                     const double vehicleDryMass )
{
    std::vector< std::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;
    terminationSettingsList.push_back( std::make_shared< PropagationTimeTerminationSettings >(
                                           simulationStartEpoch + maximumDuration ) );
    terminationSettingsList.push_back( std::make_shared< PropagationDependentVariableTerminationSettings >(
                                           std::make_shared< SingleDependentVariableSaveSettings >(
                                               altitude_dependent_variable, "Vehicle", "Moon" ), terminationAltitude, false ) );
    terminationSettingsList.push_back( std::make_shared< PropagationDependentVariableTerminationSettings >(
                                           std::make_shared< SingleDependentVariableSaveSettings >(
                                               altitude_dependent_variable, "Vehicle", "Moon" ), 0.0, true ) );
    terminationSettingsList.push_back( std::make_shared< PropagationDependentVariableTerminationSettings >(
                                           std::make_shared< SingleDependentVariableSaveSettings >(
                                               current_body_mass_dependent_variable, "Vehicle" ), vehicleDryMass, true ) );
    return std::make_shared< PropagationHybridTerminationSettings >( terminationSettingsList, true );

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
                                          altitude_dependent_variable, "Vehicle", "Moon" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          relative_speed_dependent_variable, "Vehicle", "Moon" ) );
    dependentVariablesList.push_back( std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                                          "Vehicle", flight_path_angle ) );
    return std::make_shared< DependentVariableSaveSettings >( dependentVariablesList, false );
}

void createSimulationSettings(
        NamedBodyMap& bodyMap,
        std::shared_ptr< MultiTypePropagatorSettings< double > >& propagatorSettings,
        std::shared_ptr< IntegratorSettings< > >& integratorSettings,
        const double specificImpulse )
{
    // DEFINE PROBLEM INDEPENDENT VARIABLES HERE:
    std::vector< double > decisionVariables =
    { 15629.13262285292, 21.50263026822358, -0.03344538412056863, -0.06456210720352829, 0.3943447499535977, 0.5358478897251189,
      -0.8607350478880107 };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            SIMULATION SETTINGS            /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define vehicle settings
    double vehicleMass = 4.7E3;
    double vehicleDryMass = 2.25E3;

    // Define simulation termination settings
    double maximumDuration = 86400.0;
    double terminationAltitude = 100.0E3;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT                   //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set simulation start epoch.
    double simulationStartEpoch = 0.0;

    // Create solar system bodies
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );

    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate );


    getDefaultBodySettings( bodiesToCreate );

    bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "Moon", "ECLIPJ2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationSettingsMap;

    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );

    accelerationsOfVehicle[ "Vehicle" ].push_back(
                getThrustAccelerationModelFromParameters(
                    decisionVariables, bodyMap, simulationStartEpoch, specificImpulse ) );
    accelerationSettingsMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Define propagator settings variables.
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Moon" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////   RETRIEVE DATA FOR PROPAGATION SETTINGS            ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::shared_ptr< PropagationTerminationSettings > terminationSettings = getPropagationTerminationSettings(
                simulationStartEpoch, maximumDuration, terminationAltitude, vehicleDryMass );
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave = getDependentVariableSaveSettings();
    Eigen::Vector6d systemInitialState = getInitialState( simulationStartEpoch, bodyMap );

    // Create propagator settings for mass (constant for all simulations)
    simulation_setup::SelectedMassRateModelMap massRateModelSettings;
    massRateModelSettings[ "Vehicle" ].push_back( std::make_shared< FromThrustMassModelSettings >( 1 ) );
    std::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings =
            std::make_shared< MassPropagatorSettings< double > >(
                std::vector< std::string >{ "Vehicle" }, massRateModelSettings,
                ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ), terminationSettings );

    // Define translational state propagation settings
    TranslationalPropagatorType propagatorType = cowell;
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalStatePropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationSettingsMap, bodiesToPropagate, systemInitialState,
                terminationSettings, propagatorType, dependentVariablesToSave );

    // Define full propagation settings
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList =
    { translationalStatePropagatorSettings, massPropagatorSettings };

    propagatorSettings =
            std::make_shared< MultiTypePropagatorSettings< double > >(
                propagatorSettingsList, terminationSettings, dependentVariablesToSave );


    // Create integrator settings
    integratorSettings =
            std::make_shared< IntegratorSettings< > >( rungeKutta4, simulationStartEpoch, 0.1 );
}

/*!
 *
 *   The thrust is computed based on a fixed thrust magnitude, and a variable thrust direction. The trust direction is defined
 *   on a set of 5 nodes, spread eavenly in time. At each node, a thrust angle theta is defined, which gives the angle between the
 *   -z and y angles in the ascent vehicle's vertical frame (see Mooij, 1994). Between the nodes, the thrust is linearly
 *   interpolated. If the propagation goes beyond the bounds of the nodes, the boundary value is used. The thrust profile
 *   is parameterized by the values of the vector decisionVariables
 *
 *    The entries of the vector 'decisionVariables' contains the following:
 *
 *   - Entry 0: Constant thrust magnitude
 *   - Entry 1: Constant spacing in time between nodes
 *   - Entry 2-6: Thrust angle theta, at nodes 1-5 (in order)
 *
 *   Details on the outputs written by this file can be found:
 *
 *      Benchmark data: comments for 'generateBenchmarks' function
 *      Results for integrator/propagator variations: comments under "RUN SIMULATION FOR VARIOUS SETTINGS"
 *
 */
int main( )
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );
    std::string outputPath = tudat_applications::getOutputPath( "LunarAscentDesignSpace/" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////   CREATE SIMULATION ENVIRONMENT                     ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double specificImpulse = 311.0;

    // Create simulation environment
    NamedBodyMap bodyMap;
    std::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings;
    std::shared_ptr< IntegratorSettings< > > integratorSettings;
    createSimulationSettings(
                bodyMap, propagatorSettings, integratorSettings, specificImpulse );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////   PERFORM MONTE-CARLO DESIGN SPACE EXPLORATION               //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define decision variable range
    std::vector< std::pair< double, double > > decisionVariableRange =
    {{ 5.0E3, 20.0E3 }, { 10, 100.0 }, { -0.1, 0.1 }, { -0.5, 0.5 }, { -0.7, 0.7 }, { -1.0, 1.0 }, { -1.3, 1.3 } };

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
    LunarAscentProblem lunarAscentProblem{ bodyMap, integratorSettings, propagatorSettings, decisionVariableRange, specificImpulse };

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
        lunarAscentProblem.fitness( decisionVariables );

        std::vector< double > objectives = lunarAscentProblem.getLastRunObjectives( );
        std::vector< double > constraints = lunarAscentProblem.getLastRunConstraints( );

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


//    pagmo::problem prob{ lunarAscentProblem };

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
