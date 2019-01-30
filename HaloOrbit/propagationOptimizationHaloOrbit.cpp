
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationCR3BPFullProblem.h"
#include "Tudat/Astrodynamics/Gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"
#include "Tudat/Astrodynamics/Gravitation/librationPoint.h"

#include "../applicationOutput.h"

using namespace tudat;
using namespace tudat::ephemerides;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::basic_mathematics;
using namespace tudat::orbital_element_conversions;
using namespace tudat:: propagators;

//! Create accelerations, to be used in full numerical simulation
SelectedAccelerationMap getHaloOrbitBodyMap( )
{
    SelectedAccelerationMap accelerationSettingsMap;
    accelerationSettingsMap[ "Spacecraft" ][ "Earth" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                      basic_astrodynamics::central_gravity ) );
    accelerationSettingsMap[ "Spacecraft" ][ "Sun" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                    basic_astrodynamics::central_gravity ) );

    return accelerationSettingsMap;
}


//! Create body map, to be used in simulations
NamedBodyMap getHaloOrbitBodyMap(
        const double primarySecondaryDistance,
        const double primaryGravitationalParameter,
        const double secondaryGravitationalParameter,
        const std::string& nameBodyToPropagate )
{
    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings = setupBodySettingsCR3BP(
                primarySecondaryDistance,
                "Sun", "Earth", "ECLIPJ2000", primaryGravitationalParameter, secondaryGravitationalParameter );

    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );
    bodyMap[ nameBodyToPropagate ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ nameBodyToPropagate ]->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                      std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                      < double, Eigen::Vector6d > >( ), "SSB", "ECLIPJ2000" ) );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    return bodyMap;
}

/*!
 *   This function computes the dynamics of a spacecraft, using the circularly restricted three-body problem (CR3BP), as well
 *   as a full numerical propagation, in the Earth-Sun system. The initial conditions you are provided all provide (nearly)
 *   periodic halo orbits.
 *
 *   In the code, as given to you, the dynamical model used in the numerical propagation is fully consistent with that of the
 *   CR3BP: only point mass gravity attraction by the Sun and Earth, with Sun and Earth in circular orbits around their common
 *   barycenter
 *
 *   Key outputs:
 *
 *   cr3bpStateHistory: Dynamics, in normalized, corotating coordinates, as directly computed by the CR3BP propagation.
 *   unnormalizedCr3bpStateHistory: Dynamics, in unnormalized, inertial Cartesian coordinates, as computed by the CR3BP
 *      propagation. Computed from cr3bpStateHistory in post-processing.
 *   propagatedStateHistory: Dynamics, in unnormalized, inertial Cartesian coordinates, as computed by the full
 *      numerical propagation. Directly computed by the SingleArcDynamicsSimulation
 *   normalizedPropagatedStateHistory: Dynamics, in normalized, corotating coordinates, as computed by the full
 *      numerical propagation. Computed from propagatedStateHistory in post-processing.
 *
 *   Input parameters:
 *
 *   normalizedInitialState: Initial conditions of the dynamics, given in normalized, corotating elements.
 */
int main( )
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    std::string outputPath = tudat_applications::getOutputPath( "HaloOrbit" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        ORBIT SETTINGS                 /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // DEFINE PROBLEM INDEPENDENT VARIABLES HERE:
    Eigen::Vector6d normalizedInitialState;
    normalizedInitialState <<  1.008381226326397e+00,    0.000000000000000e+00,    1.000000000000000e-04,
            0.000000000000000e+00,    9.755309237585266e-03,    0.000000000000000e+00  ;

    double primarySecondaryDistance = tudat::physical_constants::ASTRONOMICAL_UNIT;
    double primaryGravitationalParameter = simulation_setup::createGravityFieldModel(
           simulation_setup::getDefaultGravityFieldSettings( "Sun", TUDAT_NAN, TUDAT_NAN ),
            "Sun" )->getGravitationalParameter( );
    double secondaryGravitationalParameter = simulation_setup::createGravityFieldModel(
            simulation_setup::getDefaultGravityFieldSettings( "Earth", TUDAT_NAN, TUDAT_NAN ),
            "Earth" )->getGravitationalParameter( );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        OUTPUT MAPS                 ////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // State history propagated in CR3BP, unnormalized and transformed to inertial Cartesian states
    std::map< double, Eigen::Vector6d > unnormalizedCr3bpStateHistory;

    // State history propagated in CR3BP, normalized and corotating
    std::map< double, Eigen::Vector6d > cr3bpStateHistory;

    // State history numerically propagated in full dynamical model, in inertial Cartesian states
    std::map< double, Eigen::VectorXd > propagatedStateHistory;

    // State history numerically propagated in full dynamical model, converted to normalized and corotating coordinates
    std::map< double, Eigen::VectorXd > normalizedPropagatedStateHistory;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        CREATE ENVIRONMENT                 /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double initialPropagationTime = 0.0;
    double finalPropagationTime = 0.5 * tudat::physical_constants::JULIAN_YEAR;
    double integrationTimeStep = 1000.0;

    NamedBodyMap bodyMap = getHaloOrbitBodyMap(
                primarySecondaryDistance, primaryGravitationalParameter, secondaryGravitationalParameter, "Spacecraft" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////  PROPAGATE ORBIT NUMERICALLY IN CR3BP                //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double dimensionLessTimeStep =
            circular_restricted_three_body_problem::convertDimensionalTimeToDimensionlessTime(
                integrationTimeStep, primaryGravitationalParameter, secondaryGravitationalParameter,
                primarySecondaryDistance );
    double dimensionLessFinalTime = circular_restricted_three_body_problem::convertDimensionalTimeToDimensionlessTime(
                finalPropagationTime, primaryGravitationalParameter, secondaryGravitationalParameter,
                primarySecondaryDistance );
    double massParameter =
            circular_restricted_three_body_problem::computeMassParameter(
                primaryGravitationalParameter, secondaryGravitationalParameter );

    cr3bpStateHistory = performCR3BPIntegration(
                std::make_shared < numerical_integrators::IntegratorSettings < > >
                ( numerical_integrators::rungeKutta4, initialPropagationTime, dimensionLessTimeStep ),
                massParameter,  normalizedInitialState, dimensionLessFinalTime, true );

    for( auto stateIterator : cr3bpStateHistory )
    {
        unnormalizedCr3bpStateHistory[ circular_restricted_three_body_problem::convertDimensionlessTimeToDimensionalTime(
                    stateIterator.first, primaryGravitationalParameter, secondaryGravitationalParameter,
                    primarySecondaryDistance ) ] =
                circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates(
                    secondaryGravitationalParameter, primaryGravitationalParameter,
                    primarySecondaryDistance, stateIterator.second, stateIterator.first );
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////  PROPAGATE ORBIT NUMERICALLY IN FULL SYSTEM                 ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::string centralBodyOfPropagation = "SSB";

    SelectedAccelerationMap accelerationSettings  = getHaloOrbitBodyMap( );
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationSettings, { "Spacecraft" }, { centralBodyOfPropagation } );

    std::vector< std::string > centralBodies =  { centralBodyOfPropagation };
    std::vector< std::string > bodiesToPropagate = { "Spacecraft" };

    Eigen::Vector6d initialCartesianState =
            circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates(
                secondaryGravitationalParameter, primaryGravitationalParameter,
                primarySecondaryDistance, normalizedInitialState, initialPropagationTime );

    std::shared_ptr< TranslationalStatePropagatorSettings< double> > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, initialCartesianState,
              std::make_shared< PropagationTimeTerminationSettings >( finalPropagationTime, true ) );

    std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
            std::make_shared < numerical_integrators::IntegratorSettings < > >
            ( numerical_integrators::rungeKutta4, initialPropagationTime, integrationTimeStep );

    SingleArcDynamicsSimulator< > dynamicsSimulator = SingleArcDynamicsSimulator< >(
                bodyMap, integratorSettings, propagatorSettings );


    dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    for( auto stateIterator : propagatedStateHistory )
    {
        normalizedPropagatedStateHistory[ circular_restricted_three_body_problem::convertDimensionalTimeToDimensionlessTime(
                    stateIterator.first, primaryGravitationalParameter, secondaryGravitationalParameter,
                    primarySecondaryDistance ) ] =
                circular_restricted_three_body_problem::convertCartesianToCorotatingNormalizedCoordinates(
                    secondaryGravitationalParameter, primaryGravitationalParameter,
                    primarySecondaryDistance, stateIterator.second, stateIterator.first );
    }

    input_output::writeDataMapToTextFile( unnormalizedCr3bpStateHistory, "cr3bpResultUnnormalized.dat", outputPath );
    input_output::writeDataMapToTextFile( cr3bpStateHistory, "cr3bpResultNormalized.dat", outputPath );
    input_output::writeDataMapToTextFile( propagatedStateHistory, "numericalResultUnnormalized.dat", outputPath );
    input_output::writeDataMapToTextFile( normalizedPropagatedStateHistory, "numericalResultNormalized.dat", outputPath );



    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}


