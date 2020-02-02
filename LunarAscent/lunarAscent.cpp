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

LunarAscentProblem::LunarAscentProblem( const simulation_setup::NamedBodyMap bodyMap,
                                        const std::shared_ptr< IntegratorSettings< > > integratorSettings,
                                        const std::shared_ptr< PropagatorSettings< double > > propagatorSettings,
                                        const double initialTime ):
    bodyMap_(bodyMap), integratorSettings_(integratorSettings), propagatorSettings_(propagatorSettings),
    initialTime_(initialTime)
{

}

std::vector< double > LunarAscentProblem::fitness( std::vector< double >& thrustParameters ) const
{

    // Define thrust functions
    std::shared_ptr< LunarAscentThrustGuidance > thrustGuidance =
            std::make_shared< LunarAscentThrustGuidance >(
                bodyMap_.at( "Vehicle" ), initialTime_, thrustParameters );
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
    // Reset every time the fitness function is called
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                    basic_astrodynamics::central_gravity ) );

    // Should be updated during each iteration (fitness function)
    accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >(
                                                       thrustDirectionGuidanceSettings, thrustMagnitudeSettings ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    propagatorSettings_->resetIntegratedStateModels( bodyMap_ );

    dynamicsSimulator_ = std::make_shared< SingleArcDynamicsSimulator< > >( bodyMap_, integratorSettings_, propagatorSettings_ );

}
