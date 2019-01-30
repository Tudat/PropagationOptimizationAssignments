/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/Astrodynamics/TrajectoryDesign/trajectory.h>
#include <Tudat/SimulationSetup/PropagationSetup/propagationPatchedConicFullProblem.h>

#include "../applicationOutput.h"

using namespace tudat;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::gravitation;
using namespace tudat::numerical_integrators;
using namespace tudat::transfer_trajectories;

/*!
 *   This function computes a patched conic trajectory with a given set of flyby bodies, minimum periapsis distances, and
 *   transfer leg types. A Trajectory object is created that computes the required Delta V at arrival, departure, and at
 *   each flyby/DSM. Subsequently, the function fullPropagationPatchedConicsTrajectory is called, which propagates the problem
 *   numerically in the following way:
 *
 *   - The patched conic (consisting of a set of Lambert targeters), is reconstructed, and the output (vehicle state) at regular
 *     intervals is saved.
 *   - Starting at the midpoint (in time) of each transfer leg, the spacecraft state is numerically propagated using acceleration/
 *     environment settings that may be defined by the user. The propagation is done from the midpoint forwards, and backwards,
 *     to obtain the state history of the full leg.
 *
 *   In the code, as given to you, the dynamical model used in the numerical propagation is fully consistent with that of the
 *   patched conic method: only point mass gravity attraction by the Sun, with the Sun fixed at the origin.
 *
 *   Key outputs:
 *
 *   lambertTargeterResultForEachLeg: a list of the state history of the spacecraft (per leg) according to the patched conic
 *      method
 *   fullProblemResultForEachLeg: a list of the state history of the spacecraft (per leg) as produced by the numerical propagation
 *
 *   Input parameters:
 *
 *   trajectoryIndependentVariables: A vector defining the start time (first entry, in seconds since J2000) and duration of all 4
 *      legs (second, third, fourth and fifth entries respectively). For technical reasons, the final (in this case sixth) entry
 *      of this vector should always be 0.
 *   transferBodyOrder: The order of bodies (including departure and arrival bodies) for the trajectory.
 */
int main( )
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    std::string outputPath = tudat_applications::getOutputPath( "HighThrust" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        TRANSFER SETTINGS                 //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // DEFINE PROBLEM INDEPENDENT VARIABLES HERE:
    std::vector< double > trajectoryIndependentVariables =
    { 7304.5 * physical_constants::JULIAN_DAY, 300 * physical_constants::JULIAN_DAY, 300 * physical_constants::JULIAN_DAY,
      300 * physical_constants::JULIAN_DAY, 300 * physical_constants::JULIAN_DAY, 0 };

    std::vector< std::string > transferBodyOrder = { "Earth", "Venus", "Venus", "Earth", "Jupiter" };
    std::vector< TransferLegType > transferLegTypes = { mga_Departure, mga_Swingby, mga_Swingby, mga_Swingby, capture };
    std::vector< double > minimumPericenterRadii = { 1.0E4, 1.0E4, 1.0E4, 1.0E4, 1.0E4 };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        SETUP SOLAR SYSTEM BODIES            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Venus" );
    bodiesToCreate.push_back( "Jupiter" );

    NamedBodyMap bodyMap = setupBodyMapFromEphemeridesForPatchedConicsTrajectory(
                "Sun", "Spacecraft", transferBodyOrder );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE SPACECRAFT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "Spacecraft" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "Spacecraft" ]->setConstantBodyMass( 400.0 );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////             CREATE PATCHED CONIC SEMI-ANALYTICAL TRAJECTORY            ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    transfer_trajectories::Trajectory trajectory = createTransferTrajectoryObject(
                bodyMap, transferBodyOrder, "Sun", transferLegTypes, trajectoryIndependentVariables,
                minimumPericenterRadii, false, TUDAT_NAN, TUDAT_NAN, true, 1.5E9, 0.95 );

    // Retrieve total Delta V values
    double totalDeltaV;
    trajectory.calculateTrajectory( totalDeltaV );
    double captureDeltaV;
    trajectory.getCaptureDeltaV( captureDeltaV );
    double departureDeltaV;
    trajectory.getDepartureDeltaV( departureDeltaV );

    // Retrieve times, positions, and delta V at each maneuver
    std::vector < Eigen::Vector3d > positionVector;
    std::vector < double > timeVector;
    std::vector < double > deltaVVector;
    trajectory.maneuvers( positionVector, timeVector, deltaVVector );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////             NUMERICALLY PROPAGATE DYNAMICS            ////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector < basic_astrodynamics::AccelerationMap > accelerationMap =
            setupAccelerationMapPatchedConicsTrajectory(
                transferLegTypes.size( ), "Sun", "Spacecraft", bodyMap );

    std::map< int, std::map< double, Eigen::Vector6d > >& lambertTargeterResultForEachLeg;
    std::map< int, std::map< double, Eigen::Vector6d > >& fullProblemResultForEachLeg;

    std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > > differenceStateArrivalAndDeparturePerLeg =
            propagators::getDifferenceFullProblemWrtPatchedConicsTrajectory(
                bodyMap, accelerationMap, transferBodyOrder,
                "Sun", "Spacecraft", transferLegTypes, trajectoryIndependentVariables, minimumPericenterRadii,
    { TUDAT_NAN, 1.5E9 }, { TUDAT_NAN, 0.95 },
                std::make_shared < numerical_integrators::IntegratorSettings < > > (
                    numerical_integrators::rungeKutta4, 0.0, 1000.0 ), false  );

    for( auto outputIterator : differenceStateArrivalAndDeparturePerLeg )
    {
        std::cout<<outputIterator.second.first.transpose( )<<std::endl;
        std::cout<<outputIterator.second.second.transpose( )<<std::endl<<std::endl;

    }

    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
