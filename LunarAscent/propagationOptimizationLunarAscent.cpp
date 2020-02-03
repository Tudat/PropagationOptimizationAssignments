/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

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
    ascentVehicleSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 0.0;
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

//! Get the propagation termination settings for the lunar ascent.
/*!
 * \param initialTime Start time of the simulation in seconds.
 * \param maximumDuration Time in seconds, specifying the maximum time duration before which the
 * simulation should stop.
 * \param terminationAltitude Altitude in meters, specifying the maximum altitude before which the
 * simulation should stop.
 * \param vehicleDryMass Dry mass of the vehicle in kg. This is value is used to create a termination
 * condition that mandates the simulation to stop once all fuel has been used up.
 * \return Shared pointer to the PropagationTerminationSettings object.
 */
std::shared_ptr< PropagationTerminationSettings > getPropagationTerminationSettings( const double initialTime,
                                                                                     const double maximumDuration,
                                                                                     const double terminationAltitude,
                                                                                     const double vehicleDryMass )
{

    std::vector< std::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;
    terminationSettingsList.push_back( std::make_shared< PropagationTimeTerminationSettings >(
                                           initialTime + maximumDuration ) );
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
 *   This function computes the dynamics of a lunar ascent vehicle, starting at zero velocity on the Moon's surface. The only
 *   accelerations acting on the spacecraft are the Moon's point-mass gravity, and the thrust of the vehicle. Both the
 *   translational dynamics and mass of the vehicle are propagated, using a fixed specific impulse.
 *
 *   The thrust is computed based on a fixed thrust magnitude, and a variable thrust direction. The trust direction is defined
 *   on a set of 5 nodes, spread eavenly in time. At each node, a thrust angle theta is defined, which gives the angle between the
 *   -z and y angles in the ascent vehicle's vertical frame (see Mooij, 1994). Between the nodes, the thrust is linearly
 *   interpolated. If the propagation goes beyond the bounds of the nodes, the boundary value is used.
 *
 *   The propagation is terminated as soon as one of the following conditions is met:
 *
 *   - Altitude > 100 km
 *   - Altitude < 0 km
 *   - Propagation time > 3600 s
 *   - Vehicle mass < 1000 kg
 *
 *   Key outputs:
 *
 *   propagatedStateHistory Numerically propagated Cartesian state
 *   dependentVariableHistory Dependent variables saved during the state propagation of the ascent *
 *
 *   Input parameters:
 *
 *   thrustParameters: Vector contains the following:
 *
 *   - Entry 0: Constant thrust magnitude
 *   - Entry 1: Constant spacing in time between nodes
 *   - Entry 2-6: Thrust angle theta, at five nodes
 */
int main( )
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    std::string outputPath = tudat_applications::getOutputPath( "LunarAscent" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            SIMULATION SETTINGS            /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define vehicle settings
    double vehicleMass = 4.7E3;
    double vehicleDryMass = 2.25E3;
    double constantSpecificImpulse = 311.0;

    // Define simulation settings
    double initialTime = 0.0;
    double maximumDuration = 86400.0;
    double terminationAltitude = 100.0E3;



    // DEFINE PROBLEM INDEPENDENT VARIABLES HERE:
    std::vector< double > thrustParameters =
    { 15629.13262285292, 21.50263026822358, -0.03344538412056863, -0.06456210720352829, 0.3943447499535977, 0.5358478897251189,
      -0.8607350478880107 };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT                   //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create solar system bodies
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Moon" );
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
    setGlobalFrameBodyEphemerides( bodyMap, "Moon", "ECLIPJ2000" );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




    // Define propagator settings variables.
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Moon" );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationSettingsMap;

    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfVehicle[ "Vehicle" ].push_back(
                getThrustAccelerationModelFromParameters(
                        thrustParameters, bodyMap, initialTime, constantSpecificImpulse ) );
    accelerationSettingsMap[ "Vehicle" ] = accelerationsOfVehicle;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Define dependent variables
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          altitude_dependent_variable, "Vehicle", "Moon" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          relative_speed_dependent_variable, "Vehicle", "Moon" ) );
    dependentVariablesList.push_back( std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                                          "Vehicle", flight_path_angle ) );

    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

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

    std::shared_ptr< PropagationTerminationSettings > terminationSettings = getPropagationTerminationSettings( initialTime,
                                                                                                               maximumDuration,
                                                                                                               terminationAltitude,
                                                                                                               vehicleDryMass );

    // Loop over settings
    const unsigned int numberOfPropagators = 7;
    const unsigned int numberOfIntegrators = 5;
    int numberOfIntegratorStepSizeSettings = 4;
    for( unsigned int i = 0; i < numberOfPropagators; i++ )
    {

        for( unsigned int j = 0; j < numberOfIntegrators; j++ )
        {

            for( int k = 0; k < numberOfIntegratorStepSizeSettings; k++ )
            {
                // Define propagator type
                TranslationalPropagatorType propagatorType = cowell;

                // Retrieve initial state
                Eigen::Vector6d systemInitialState = getInitialState( initialTime, bodyMap );

                // Define translational state propagation settings
                std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalStatePropagatorSettings =
                        std::make_shared< TranslationalStatePropagatorSettings< double > >(
                            centralBodies, accelerationSettingsMap, bodiesToPropagate, systemInitialState,
                            terminationSettings, propagatorType );

                // Define mass propagation settings

                simulation_setup::SelectedMassRateModelMap massRateModelSettings;
                massRateModelSettings[ "Vehicle" ].push_back( std::make_shared< FromThrustMassModelSettings >( 1 ) );

                std::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings =
                        std::make_shared< MassPropagatorSettings< double > >(
                            std::vector< std::string >{ "Vehicle" }, massRateModelSettings,
                            ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ), terminationSettings );

                // Define full propagation settings
                std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector =
                { translationalStatePropagatorSettings, massPropagatorSettings };

                std::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings =
                        std::make_shared< MultiTypePropagatorSettings< double > >(
                            propagatorSettingsVector, terminationSettings, dependentVariablesToSave );

                // Define integration settings
                std::shared_ptr< IntegratorSettings< > > integratorSettings =
                        std::make_shared< IntegratorSettings< > >( rungeKutta4, initialTime, 1.0 );

                ////////////////////////////////////~///////////////////////////////////////////////////////////////////////////////////
                ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                LunarAscentProblem prob{ bodyMap, integratorSettings, propagatorSettings, initialTime, constantSpecificImpulse };

                prob.fitness( thrustParameters );

                std::map< double, Eigen::VectorXd > propagatedStateHistory = prob.getLastRunPropagatedStateHistory();
                std::map< double, Eigen::VectorXd > dependentVariableHistory = prob.getLastRunDependentVariableHistory();

                input_output::writeDataMapToTextFile( propagatedStateHistory, "stateHistory.dat", outputPath );
                input_output::writeDataMapToTextFile( dependentVariableHistory, "dependentVariables.dat", outputPath );
                input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                                     thrustParameters ), "thrustParameters.dat", 16, outputPath );
            }

        }

    }



    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
