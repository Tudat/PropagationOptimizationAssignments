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

namespace tudat_applications
{
namespace PropagationOptimization2020
{

double getTrajectoryTimeOfFlight( const std::vector< double >& decisionVariables)
{
    return decisionVariables.at( 1 ) * physical_constants::JULIAN_DAY;
}

double getTrajectoryInitialTime( const std::vector< double >& decisionVariables, const double bufferTime  )
{
    return decisionVariables.at( 0 ) * physical_constants::JULIAN_DAY + bufferTime;
}

double getTrajectoryFinalTime( const std::vector< double >& decisionVariables, const double bufferTime )
{
    return ( decisionVariables.at( 0 ) + decisionVariables.at( 1 )  )* physical_constants::JULIAN_DAY - bufferTime;
}


//! Get the propagation termination settings for the lunar ascent.
/*!
 * \param simulationStartEpoch Start time of the simulation in seconds.
 * \param maximumDuration Time in seconds, specifying the maximum time duration before which the
 * simulation should stop.
 * \param terminationAltitude Altitude in meters, specifying the maximum altitude before which the
 * simulation should stop.
 * \param vehicleDryMass Dry mass of the vehicle in kg. This is value is used to create a termination
 * condition that mandates the simulation to stop once all fuel has been used up.
 * \return Shared pointer to the PropagationTerminationSettings object.
 */
std::shared_ptr< PropagationTerminationSettings > getPropagationTerminationSettings(
        const std::vector< double >& decisionVariables,
        const double targetDistance,
        const double trajectoryTimeBuffer )
{

    std::vector< std::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;
    terminationSettingsList.push_back( std::make_shared< PropagationTimeTerminationSettings >(
                                           getTrajectoryFinalTime( decisionVariables, trajectoryTimeBuffer ) ) );
    terminationSettingsList.push_back( std::make_shared< PropagationDependentVariableTerminationSettings >(
                                           std::make_shared< SingleDependentVariableSaveSettings >(
                                               relative_distance_dependent_variable, "Vehicle", "Mars" ), targetDistance, true ) );
    return std::make_shared< PropagationHybridTerminationSettings >( terminationSettingsList, true );

}

void getRadialVelocityShapingFunctions(
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& radialVelocityFunctionComponents,
        Eigen::VectorXd& freeCoefficientsRadialVelocityFunction,
        const std::vector< double >& decisionVariables, const double frequency, const double scaleFactor,
        const double timeOfFlight, const int numberOfRevolutions  )
{
    // Retrieve default methods (lowest-order in Gondelach and Noomen, 2015)
    shape_based_methods::getRecommendedRadialVelocityBaseFunctions(
                radialVelocityFunctionComponents, freeCoefficientsRadialVelocityFunction, timeOfFlight );

    // Add degrees of freedom (highest-order in Gondelach and Noomen, 2015)
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > fourthRadialVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                1.0, 0.5 * frequency, scaleFactor );
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > fifthRadialVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                1.0, 0.5 * frequency, scaleFactor );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerSine, fourthRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerCosine, fifthRadialVelocityBaseFunctionSettings ) );

    // Set free parameters
    freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 2 );
    freeCoefficientsRadialVelocityFunction( 0 ) = decisionVariables.at( 3 );
    freeCoefficientsRadialVelocityFunction( 1 ) = decisionVariables.at( 4 );
}

void getNormalVelocityShapingFunctions(
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& normalVelocityFunctionComponents,
        Eigen::VectorXd& freeCoefficientsNormalVelocityFunction,
        const std::vector< double >& decisionVariables, const double frequency, const double scaleFactor,
        const double timeOfFlight, const int numberOfRevolutions  )
{
    // Retrieve default methods (lowest-order in Gondelach and Noomen, 2015)
    shape_based_methods::getRecommendedNormalAxialBaseFunctions( normalVelocityFunctionComponents, freeCoefficientsNormalVelocityFunction, timeOfFlight );

    // Add degrees of freedom (highest-order in Gondelach and Noomen, 2015)
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > fourthNormalVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                1.0, 0.5 * frequency, scaleFactor );
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > fifthNormalVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                1.0, 0.5 * frequency, scaleFactor );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerSine, fourthNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerCosine, fifthNormalVelocityBaseFunctionSettings ) );

    // Set free parameters
    freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 2 );
    freeCoefficientsNormalVelocityFunction( 0 ) = decisionVariables.at( 5 );
    freeCoefficientsNormalVelocityFunction( 1 ) = decisionVariables.at( 6 );
}

void getAxialVelocityShapingFunctions(
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& axialVelocityFunctionComponents,
        Eigen::VectorXd& freeCoefficientsAxialVelocityFunction,
        const std::vector< double >& decisionVariables, const double frequency, const double scaleFactor,
        const double timeOfFlight, const int numberOfRevolutions )
{
    // Retrieve default methods (lowest-order in Gondelach and Noomen, 2015)
    shape_based_methods::getRecommendedAxialVelocityBaseFunctions(
                axialVelocityFunctionComponents, freeCoefficientsAxialVelocityFunction,
                timeOfFlight, numberOfRevolutions );

    // Add degrees of freedom (highest-order in Gondelach and Noomen, 2015)
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > fourthAxialVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >
            ( 3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > fifthAxialVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerCosine, fourthAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerSine, fifthAxialVelocityBaseFunctionSettings ) );

    // Set free parameters
    freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 2 );
    freeCoefficientsAxialVelocityFunction( 0 ) = decisionVariables.at( 7 );
    freeCoefficientsAxialVelocityFunction( 1 ) =  decisionVariables.at( 8 );
}



//! Function that creates the object that computes the semi-analytical hodographic shaping method
std::shared_ptr< HodographicShaping > createHodographicShapingObject(
        const std::vector< double >& decisionVariables,
        const simulation_setup::NamedBodyMap bodyMap )
{
    // Retrieve basic properties from trajectory parameters
    double initialTime = getTrajectoryInitialTime( decisionVariables );
    double timeOfFlight = getTrajectoryTimeOfFlight( decisionVariables );
    double finalTime = getTrajectoryFinalTime( decisionVariables );
    int numberOfRevolutions = decisionVariables.at( 2 );

    // Compute relevant frequency and scale factor for shaping functions
    double frequency = 2.0 * mathematical_constants::PI / timeOfFlight;
    double scaleFactor = 1.0 / timeOfFlight;

    // Retrieve settings for base functions in radial, normal and axial directions
    Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 0 );
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );

    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;

    getRadialVelocityShapingFunctions(
                radialVelocityFunctionComponents, freeCoefficientsRadialVelocityFunction, decisionVariables,
                frequency, scaleFactor, timeOfFlight, numberOfRevolutions );
    getNormalVelocityShapingFunctions(
                normalVelocityFunctionComponents, freeCoefficientsNormalVelocityFunction, decisionVariables,
                frequency, scaleFactor, timeOfFlight, numberOfRevolutions );
    getAxialVelocityShapingFunctions(
                axialVelocityFunctionComponents, freeCoefficientsAxialVelocityFunction, decisionVariables,
                frequency, scaleFactor, timeOfFlight, numberOfRevolutions );

    // Retrieve boundary conditions and central body gravitational parameter
    Eigen::Vector6d initialState = bodyMap.at( "Earth" )->getStateInBaseFrameFromEphemeris( initialTime );
    Eigen::Vector6d finalState = bodyMap.at( "Mars" )->getStateInBaseFrameFromEphemeris( finalTime );
    double centralBodyGravitationalParameter = bodyMap.at( "Sun" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Create and return shape-based method
    return std::make_shared< HodographicShaping >(
                initialState, finalState, timeOfFlight, centralBodyGravitationalParameter, numberOfRevolutions,
                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction,
                freeCoefficientsNormalVelocityFunction,
                freeCoefficientsAxialVelocityFunction );
}

//! Function that generates thrust acceleration model from thrust parameters
std::shared_ptr< ThrustAccelerationSettings > getThrustAccelerationSettingsFromParameters(
        const std::vector< double >& decisionVariables,
        const simulation_setup::NamedBodyMap bodyMap )
{
    return createHodographicShapingObject( decisionVariables, bodyMap )->getLowThrustAccelerationSettings(
                bodyMap, "Vehicle", nullptr, nullptr, getTrajectoryInitialTime( decisionVariables ) );
}

//! Function to retrieve the semi-analytical calculation of the state along the shape-based trajectory at a given time
Eigen::Vector6d getHodographicLowThrustStateAtEpoch(
        const std::vector< double >& decisionVariables,
        const simulation_setup::NamedBodyMap bodyMap,
        const double evaluationTime )
{
    std::map< double, Eigen::Vector6d > stateMap;
    std::vector< double > epochs = { evaluationTime - getTrajectoryInitialTime( decisionVariables, 0.0 ) };
    createHodographicShapingObject( decisionVariables, bodyMap )->getTrajectory( epochs, stateMap );
    return stateMap.begin( )->second;

}

}

}

using namespace tudat_applications::PropagationOptimization2020;

//! Constructor for the problem class
LowThrustProblem::LowThrustProblem( const simulation_setup::NamedBodyMap bodyMap,
                                    const std::shared_ptr< IntegratorSettings< > > integratorSettings,
                                    const std::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings,
                                    const std::vector< std::pair< double, double > >& decisionVariableRange,
                                    double specificImpulse,
                                    double minimumMarsDistance,
                                    double timeBuffer ):
    bodyMap_(bodyMap), integratorSettings_(integratorSettings), propagatorSettings_(propagatorSettings),
    specificImpulse_( specificImpulse ), minimumMarsDistance_( minimumMarsDistance ), timeBuffer_( timeBuffer )
{
    translationalStatePropagatorSettings_ =
            std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< double > >(
                propagatorSettings_->propagatorSettingsMap_.at( translational_state ).at( 0 ) );

    std::vector< double > boxBoundMinima, boxBoundMaxima;
    for( unsigned int i = 0; i < decisionVariableRange.size( ); i++ )
    {
        boxBoundMinima.push_back( decisionVariableRange.at( i ).first );
        boxBoundMaxima.push_back( decisionVariableRange.at( i ).second );
    }
    boxBounds_ = std::make_pair( boxBoundMinima, boxBoundMaxima );
}

//! Function to compute propagate the dynamics of the vehicle defined by the decisionVariables
std::vector< double > LowThrustProblem::fitness( const std::vector< double >& decisionVariables ) const
{
    // Create trajectory shape object
    hodographicShaping_ = createHodographicShapingObject(
                decisionVariables, bodyMap_ );

    double currentDeltaV = hodographicShaping_->computeDeltaV( );
    double deltaVLimitValue = 50.0E3;

    bool performPropagation = ( currentDeltaV < deltaVLimitValue );

    std::cout<<"Calc: "<<currentDeltaV<<std::endl;
    if( performPropagation )
    {
        double initialPropagationTime = getTrajectoryInitialTime( decisionVariables, timeBuffer_ );
        integratorSettings_->initialTime_ = initialPropagationTime;

        // Extract existing acceleration settings, and clear existing self-exerted accelerations of vehicle
        simulation_setup::SelectedAccelerationMap accelerationSettings =
                translationalStatePropagatorSettings_->getAccelerationSettingsMap( );
        accelerationSettings[ "Vehicle" ][ "Vehicle" ].clear( );

        // Retrieve new acceleration model for thrust and set in list of settings
        std::shared_ptr< AccelerationSettings > newThrustSettings = hodographicShaping_->getLowThrustAccelerationSettings(
                    bodyMap_, "Vehicle", [=](const double){return specificImpulse_;}, integratorSettings_,
        getTrajectoryInitialTime( decisionVariables ) );
        accelerationSettings[ "Vehicle" ][ "Vehicle" ].push_back( newThrustSettings );

        // Update translational propagatot settings
        translationalStatePropagatorSettings_->resetAccelerationModelsMap(
                    accelerationSettings, bodyMap_ );
        Eigen::Vector6d systemInitialState = getHodographicLowThrustStateAtEpoch(
                    decisionVariables, bodyMap_, initialPropagationTime );
        translationalStatePropagatorSettings_->resetInitialStates( systemInitialState );

        // Update full propagator settings
        propagatorSettings_->resetIntegratedStateModels( bodyMap_ );
        propagatorSettings_->resetInitialStates(
                    createCombinedInitialState< double >( propagatorSettings_->propagatorSettingsMap_ ).segment( 0, 7 ) );
        propagatorSettings_->resetTerminationSettings(
                    getPropagationTerminationSettings(
                        decisionVariables, minimumMarsDistance_, 0.0 ) );

        dynamicsSimulator_ = std::make_shared< SingleArcDynamicsSimulator< > >( bodyMap_, integratorSettings_, propagatorSettings_ );

    }

    computeObjectivesAndConstraints( decisionVariables, performPropagation );
    return objectives_;

}

void LowThrustProblem::computeObjectivesAndConstraints( const std::vector< double >& decisionVariables,
                                                        const bool performPropagation ) const
{
    std::map< double, Eigen::VectorXd > stateHistory;
    std::map< double, Eigen::VectorXd > dependentVariableHistory;

    if( performPropagation )
    {
        stateHistory = dynamicsSimulator_->getEquationsOfMotionNumericalSolution( );
        dependentVariableHistory = dynamicsSimulator_->getDependentVariableHistory( );
    }

    double currentDeltaV = hodographicShaping_->computeDeltaV( );

    constraints_; // =
    objectives_ = { currentDeltaV }; // =
}
