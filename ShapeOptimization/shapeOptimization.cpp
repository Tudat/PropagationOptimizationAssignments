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

namespace tudat_applications
{
namespace PropagationOptimization2020
{
/*!
 *  Function that creates an aerodynamic database for a capsule, based on a set of shape parameters
 *  The Capsule shape consists of four separate geometrical components: a sphere segment for the nose, a torus segment for the
 *  shoulder/edge, a conical frustum for the rear body, and a sphere segment for the rear cap (see Dirkx and Mooij, 2016).
 *  The code used in this function discretizes these surfaces into a structured mesh of quadrilateral panels. The parameters
 *  numberOfPoints and numberOfLines define the number of discretization points (for each part) in both independent directions
 *  (lengthwise and circumferential). The list selectedMethods defines the type of aerodynamic analysis method that is used.
 */
std::shared_ptr< HypersonicLocalInclinationAnalysis > getCapsuleCoefficientInterface(
        const std::shared_ptr< geometric_shapes::Capsule > capsule,
        const std::string directory,
        const std::string filePrefix )
{

    // Define settings for surface discretization of capsule
    std::vector< int > numberOfLines;
    std::vector< int > numberOfPoints;
    numberOfLines.resize( 4 );
    numberOfPoints.resize( 4 );
    numberOfLines[ 0 ] = 31;
    numberOfPoints[ 0 ] = 31;
    numberOfLines[ 1 ] = 31;
    numberOfPoints[ 1 ] = 31;
    numberOfLines[ 2 ] = 31;
    numberOfPoints[ 2 ] = 31;
    numberOfLines[ 3 ] = 11;
    numberOfPoints[ 3 ] = 11;

    // DO NOT CHANGE THESE (setting to true will turn parts of the vehicle 'inside out')
    std::vector< bool > invertOrders;
    invertOrders.resize( 4 );
    invertOrders[ 0 ] = 0;
    invertOrders[ 1 ] = 0;
    invertOrders[ 2 ] = 0;
    invertOrders[ 3 ] = 0;

    // Define moment reference point
    Eigen::Vector3d momentReference;
    momentReference( 0 ) = -0.6624;
    momentReference( 1 ) = 0.0;
    momentReference( 2 ) = 0.1369;

    // Define independent variable values
    std::vector< std::vector< double > > independentVariableDataPoints;
    independentVariableDataPoints.resize( 3 );
    independentVariableDataPoints[ 0 ] = getDefaultHypersonicLocalInclinationMachPoints( "Full" );
    std::vector< double > angleOfAttackPoints;
    angleOfAttackPoints.resize( 15 );
    for ( int i = 0; i < 15; i++ )
    {
        angleOfAttackPoints[ i ] = static_cast< double >( i - 6 ) * 5.0 * PI / 180.0;
    }
    independentVariableDataPoints[ 1 ] = angleOfAttackPoints;
    independentVariableDataPoints[ 2 ] =
            getDefaultHypersonicLocalInclinationAngleOfSideslipPoints( );

    // Define local inclination methods to use
    std::vector< std::vector< int > > selectedMethods;
    selectedMethods.resize( 2 );
    selectedMethods[ 0 ].resize( 4 );
    selectedMethods[ 1 ].resize( 4 );

    selectedMethods[ 0 ][ 0 ] = 0;
    selectedMethods[ 0 ][ 1 ] = 0;
    selectedMethods[ 0 ][ 2 ] = 0;
    selectedMethods[ 0 ][ 3 ] = 0;
    selectedMethods[ 1 ][ 0 ] = 0;
    selectedMethods[ 1 ][ 1 ] = 0;
    selectedMethods[ 1 ][ 2 ] = 0;
    selectedMethods[ 1 ][ 3 ] = 0;

    // Create aerodynamic database
    std::shared_ptr< HypersonicLocalInclinationAnalysis > hypersonicLocalInclinationAnalysis =
            std::make_shared< HypersonicLocalInclinationAnalysis >(
                independentVariableDataPoints, capsule, numberOfLines, numberOfPoints,
                invertOrders, selectedMethods, PI * std::pow( capsule->getMiddleRadius( ), 2.0 ),
                capsule->getMiddleRadius( ), momentReference, false );

    // Save vehicle mesh to a file
    aerodynamics::saveVehicleMeshToFile(
                hypersonicLocalInclinationAnalysis, directory, filePrefix );

    // Create analysis object and capsule database.
    return  hypersonicLocalInclinationAnalysis;
}

//! Function to set vehicle properties related to vehicle shape (mass, aerodynamic coefficients)
void setVehicleShapeParameters(
        std::vector< double > shapeParameters,
        const NamedBodyMap& bodyMap,
        const double vehicleDensity )
{
    // Apply shape constraint
    double limitLength = ( shapeParameters[ 1 ] - shapeParameters[ 4 ] * ( 1.0 - std::cos( shapeParameters[ 3 ] ) ) ) /
            std::tan( -shapeParameters[ 3 ] );
    if( shapeParameters[ 2 ] >= limitLength - 0.01 )
    {
        shapeParameters[ 2 ] = limitLength - 0.01;
    }

    // Create capsule.
    std::shared_ptr< geometric_shapes::Capsule > capsule
            = std::make_shared< geometric_shapes::Capsule >(
                shapeParameters[ 0 ], shapeParameters[ 1 ], shapeParameters[ 2 ],
            shapeParameters[ 3 ], shapeParameters[ 4 ] );

    // Vehicle properties
    bodyMap.at( "Capsule" )->setConstantBodyMass(
                capsule->getVolume( ) * vehicleDensity );

    // Create vehicle aerodynamic coefficients
    bodyMap.at( "Capsule" )->setAerodynamicCoefficientInterface(
                getCapsuleCoefficientInterface(
                    capsule, tudat_applications::getOutputPath( "ShapeOptimization" ), "output_" ) );


}

//! Function to add the Capsule, and associated shape properties, to the body map.
void addCapsuleToBodyMap( NamedBodyMap& bodyMap,
                          std::vector< double >& shapeParameters,
                          const double vehicleDensity  )
{
    // Create vehicle objects.
    bodyMap[ "Capsule" ] = std::make_shared< simulation_setup::Body >( );
    setVehicleShapeParameters( shapeParameters, bodyMap, vehicleDensity );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "Earth", "J2000" );

}

}

}


//! Function to compute propagate the dynamics of the capsule defined by given shapeParameters, and compute its fitness (undefined)
std::vector< double > ShapeOptimizationProblem::fitness( std::vector< double >& shapeParameters) const
{
    // Recreate capsule body with new shape parameters
    bodyMap_.at( "Capsule" ).reset( );
    addCapsuleToBodyMap( bodyMap_, shapeParameters, vehicleDensity_ );

    // Reset propagation and guidance models
    propagatorSettings_->resetIntegratedStateModels( bodyMap_ );
    std::shared_ptr< CapsuleAerodynamicGuidance > capsuleGuidance =
            std::make_shared< CapsuleAerodynamicGuidance >( bodyMap_, shapeParameters.at( 5 ) );
    setGuidanceAnglesFunctions( capsuleGuidance, bodyMap_.at( "Capsule" ) );

    // Propagate dynamics
    dynamicsSimulator_ = std::make_shared< SingleArcDynamicsSimulator< > >( bodyMap_, integratorSettings_, propagatorSettings_ );

    // Return fitness; for now just return {0.0}
    return {0.0};

}
