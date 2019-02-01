/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/Astrodynamics/Aerodynamics/hypersonicLocalInclinationAnalysis.h>
#include <Tudat/Mathematics/GeometricShapes/capsule.h>
#include <Tudat/Mathematics/Statistics/randomVariableGenerator.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "../applicationOutput.h"

using namespace tudat::ephemerides;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::propagators;
using namespace tudat::aerodynamics;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace tudat::mathematical_constants;
using namespace tudat::reference_frames;
using namespace tudat;

/*!
 *  Class to compute the thrust direction and magnitude for the lunar ascent vehicle. The current inputs set a
 *  constant thrust magnitude, and a thrust direction in y-(-z) plane of the vertical frame define linearly in time using
 *  equispaced nodes. These settings are to be modified for the assignment.
 */
class LunarAscentThrustGuidance
{
public:

    //! Contructor
    /*!
     * Contructor
     * \param vehicleBody Body object for the ascent vehicle
     * \param initialTime Start time of the propagatiin
     * \param parameterVector Vector of independent variables to be used for thrust parameterization
     * (see main function documentation)
     */
    LunarAscentThrustGuidance(
            const std::shared_ptr< Body > vehicleBody,
            const double initialTime,
            const std::vector< double > parameterVector ):
        vehicleBody_( vehicleBody ),
        parameterVector_( parameterVector )
    {
        // Retrieve parameters of thrust profile
        thrustMagnitude_ = parameterVector_.at( 0 );
        timeInterval_ = parameterVector_.at( 1 );

        // Create interpolator for thrust angle
        double currentTime = initialTime;
        for( unsigned int i = 0; i < parameterVector_.size( ) - 2; i++ )
        {
            thrustAngleMap_[ currentTime ] = parameterVector_.at( i + 2 );
            currentTime += timeInterval_;
        }
        thrustAngleInterpolator_ = createOneDimensionalInterpolator(
                    thrustAngleMap_, std::make_shared< InterpolatorSettings >(
                        linear_interpolator, huntingAlgorithm, false, use_boundary_value ) );
    }

    //! Function that computes the inertial thrust direction for each state derivative function evaluation
    Eigen::Vector3d getCurrentThrustDirection( const double currentTime )
    {
        if(  vehicleFlightConditions_ == nullptr )
        {
            vehicleFlightConditions_ = vehicleBody_->getFlightConditions( );
        }

        // Retrieve thrust angle
        double currentThrustAngle = thrustAngleInterpolator_->interpolate( currentTime );

        // Set thrust in V-frame
        Eigen::Vector3d thrustDirectionInVerticalFrame =
                ( Eigen::Vector3d( ) << 0.0, std::sin( currentThrustAngle ), -std::cos( currentThrustAngle ) ).finished( );

        // Retrieve rotation from V-frame to inertial frame
        Eigen::Quaterniond verticalToInertialFrame =
                vehicleFlightConditions_->getAerodynamicAngleCalculator( )->getRotationQuaternionBetweenFrames(
                    vertical_frame, inertial_frame );

        // Return thrust direction
        return verticalToInertialFrame * thrustDirectionInVerticalFrame;

    }

    //! Function that computes the thrust magnitude for each state derivative function evaluation
    double getCurrentThrustMagnitude( const double currentTime )
    {
        return thrustMagnitude_;
    }

private:

    std::shared_ptr< Body > vehicleBody_;

    std::vector< double > parameterVector_;

    std::map< double, double > thrustAngleMap_;

    std::shared_ptr< OneDimensionalInterpolator< double, double > > thrustAngleInterpolator_;

    std::shared_ptr< FlightConditions > vehicleFlightConditions_;

    double timeInterval_;

    double thrustMagnitude_;

};

/*!
 *   This function computes the dynamics of a lunar ascent vehicle, starting at zero velocity on the Moon's surface. The only
 *   accelerations acting on the spacecraft are the Moon's point-mass gravity, and the thrust of the vehicle. Both the
 *   translational dynamics and mass of the vehicle are propagated, using a specific impulse Isp of 250 s.
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

    // DEFINE PROBLEM INDEPENDENT VARIABLES HERE:
    std::vector< std::pair< double, double > > thrustParametersLimits =
        {{ 5.0E3, 20.0E3 }, { 10, 100.0 }, { -0.1, 0.1 }, { -0.5, 0.5 }, { -0.7, 0.7 }, { -1.0, 1.0 }, { -1.3, 1.3 } };
    std::vector< std::function< double( ) > > parameterMonteCarloFunctions;
    for( int i = 0; i < thrustParametersLimits.size( ); i++ )
    {
        parameterMonteCarloFunctions.push_back(
                    statistics::createBoostContinuousRandomVariableGeneratorFunction(
                        statistics::uniform_boost_distribution,
        { thrustParametersLimits.at( i ).first, thrustParametersLimits.at( i ).second }, i ) );
    }


    for( int i = 0; i < 1000; i++ )
    {
        std::cout<<i<<std::endl;
        std::vector< double > thrustParameters;
        for( int j = 0; j < parameterMonteCarloFunctions.size( ); j++ )
        {
            thrustParameters.push_back( parameterMonteCarloFunctions.at( j )( ) );
        }
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Define thrust functions
        std::shared_ptr< LunarAscentThrustGuidance > thrustGuidance =
                std::make_shared< LunarAscentThrustGuidance >(
                    bodyMap.at( "Vehicle" ), initialTime, thrustParameters );
        std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction =
                std::bind( &LunarAscentThrustGuidance::getCurrentThrustDirection, thrustGuidance, std::placeholders::_1 );
        std::function< double( const double ) > thrustMagnitudeFunction =
                std::bind( &LunarAscentThrustGuidance::getCurrentThrustMagnitude, thrustGuidance, std::placeholders::_1 );

        std::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings =
                std::make_shared< CustomThrustDirectionSettings >( thrustDirectionFunction );
        std::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings =
                std::make_shared< FromFunctionThrustMagnitudeSettings >(
                    thrustMagnitudeFunction, [ = ]( const double ){ return constantSpecificImpulse; } );

        // Define acceleration settings
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
        accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                        basic_astrodynamics::central_gravity ) );
        accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >(
                                                           thrustDirectionGuidanceSettings, thrustMagnitudeSettings ) );
        accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;


        // Define propagator settings variables.
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;
        bodiesToPropagate.push_back( "Vehicle" );
        centralBodies.push_back( "Moon" );

        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // Convert the state to the global (inertial) frame.
        std::shared_ptr< ephemerides::RotationalEphemeris > moonRotationalEphemeris =
                bodyMap.at( "Moon" )->getRotationalEphemeris( );
        Eigen::VectorXd systemInitialState = transformStateToGlobalFrame(
                    bodyFixedSystemInitialState, initialTime, moonRotationalEphemeris );

        // Create termination settings.
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
        std::shared_ptr< PropagationTerminationSettings > terminationSettings = std::make_shared<
                PropagationHybridTerminationSettings >( terminationSettingsList, true );

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

        // Define translational state propagation settings
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalStatePropagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                    terminationSettings );

        // Define mass propagation settings
        std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
        massRateModels[ "Vehicle" ] = (
                    createMassRateModel( "Vehicle", std::make_shared< FromThrustMassModelSettings >( 1 ),
                                         bodyMap, accelerationModelMap ) );
        std::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings =
                std::make_shared< MassPropagatorSettings< double > >(
                    std::vector< std::string >{ "Vehicle" }, massRateModels,
                    ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ), terminationSettings );

        // Define full propagation settings
        std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector =
        { translationalStatePropagatorSettings, massPropagatorSettings };
        std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
                std::make_shared< MultiTypePropagatorSettings< double > >(
                    propagatorSettingsVector, terminationSettings, dependentVariablesToSave );

        // Define integration settings
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >( rungeKutta4, initialTime, 1.0 );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings );

        std::map< double, Eigen::VectorXd > propagatedStateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );

        input_output::writeDataMapToTextFile( propagatedStateHistory, "stateHistory" + std::to_string( i ) + ".dat", outputPath );
        input_output::writeDataMapToTextFile( dependentVariableHistory, "dependentVariables" + std::to_string( i ) + ".dat", outputPath );
        input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                             thrustParameters ), "thrustParameters" + std::to_string( i ) + ".dat", 16, outputPath );

    }
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
