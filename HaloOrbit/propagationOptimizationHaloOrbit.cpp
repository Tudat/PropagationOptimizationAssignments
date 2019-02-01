
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
SelectedAccelerationMap getHaloOrbitAccelerationsMap( )
{
    SelectedAccelerationMap accelerationSettingsMap;
    accelerationSettingsMap[ "Spacecraft" ][ "Earth" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                      basic_astrodynamics::central_gravity ) );
    accelerationSettingsMap[ "Spacecraft" ][ "Sun" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                    basic_astrodynamics::central_gravity ) );
    accelerationSettingsMap[ "Spacecraft" ][ "Mars" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                     basic_astrodynamics::central_gravity ) );
    accelerationSettingsMap[ "Spacecraft" ][ "Jupiter" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
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
    bodySettings[ "Mars" ] = simulation_setup::getDefaultSingleBodySettings( "Mars", TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );
    bodySettings[ "Jupiter" ] = simulation_setup::getDefaultSingleBodySettings( "Jupiter", TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );

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

    Eigen::MatrixXd initialConditionsList = input_output::readMatrixFromFile(
                "/home/dominic/Software/tudatBundleTest/tudatBundle/tudatApplications/PropagationOptimizationAssignments/HaloOrbit/L2_2_initial_conditions_ES.txt" );

    std::cout<<"Size: "<<initialConditionsList.rows( )<<" "<<initialConditionsList.cols( )<<std::endl;

    for( int i = 0; i < initialConditionsList.rows( ); i++ )
    {
        std::cout<<i<<std::endl;
        // DEFINE PROBLEM INDEPENDENT VARIABLES HERE:
        Eigen::Vector6d normalizedInitialState =
                initialConditionsList.block( i, 1, 1, 6 ).transpose( );

        double primarySecondaryDistance = tudat::physical_constants::ASTRONOMICAL_UNIT;
        double primaryGravitationalParameter = simulation_setup::createGravityFieldModel(
                    simulation_setup::getDefaultGravityFieldSettings( "Sun", TUDAT_NAN, TUDAT_NAN ),
                    "Sun" )->getGravitationalParameter( );
        double secondaryGravitationalParameter = simulation_setup::createGravityFieldModel(
                    simulation_setup::getDefaultGravityFieldSettings( "Earth", TUDAT_NAN, TUDAT_NAN ),
                    "Earth" )->getGravitationalParameter( );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////        CREATE ENVIRONMENT                 /////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        double initialTotalPropagationTime = 0.0 * tudat::physical_constants::JULIAN_YEAR;
        double finalTotalPropagationTime = 2.0 * tudat::physical_constants::JULIAN_YEAR;
        double integrationTimeStep = 1.0E3;

        int numberOfArcs = 12;
        double arcDuration = ( finalTotalPropagationTime - initialTotalPropagationTime ) /
                static_cast< double >( numberOfArcs );

        NamedBodyMap bodyMap = getHaloOrbitBodyMap(
                    primarySecondaryDistance, primaryGravitationalParameter, secondaryGravitationalParameter, "Spacecraft" );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////  PROPAGATE ORBIT NUMERICALLY IN CR3BP                //////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        Eigen::VectorXd currentNormalizedInitialState = normalizedInitialState;

        for( int j = 0; j < numberOfArcs; j++ )
        {

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

            double initialPropagationTime = initialTotalPropagationTime +
                    static_cast< double >( j ) * arcDuration;
            double finalPropagationTime = initialTotalPropagationTime +
                    static_cast< double >( j + 1 ) * arcDuration;

            std::cout<<j<<" "<<initialPropagationTime<<" "<<finalPropagationTime<<std::endl;

            double dimensionLessTimeStep =
                    circular_restricted_three_body_problem::convertDimensionalTimeToDimensionlessTime(
                        integrationTimeStep, primaryGravitationalParameter, secondaryGravitationalParameter,
                        primarySecondaryDistance );
            double dimensionLessInitialTime = circular_restricted_three_body_problem::convertDimensionalTimeToDimensionlessTime(
                        initialPropagationTime, primaryGravitationalParameter, secondaryGravitationalParameter,
                        primarySecondaryDistance );
            double dimensionLessFinalTime = circular_restricted_three_body_problem::convertDimensionalTimeToDimensionlessTime(
                        finalPropagationTime, primaryGravitationalParameter, secondaryGravitationalParameter,
                        primarySecondaryDistance );
            double massParameter =
                    circular_restricted_three_body_problem::computeMassParameter(
                        primaryGravitationalParameter, secondaryGravitationalParameter );

            cr3bpStateHistory = performCR3BPIntegration(
                        std::make_shared < numerical_integrators::IntegratorSettings < > >
                        ( numerical_integrators::rungeKutta4, dimensionLessInitialTime, dimensionLessTimeStep ),
                        massParameter, currentNormalizedInitialState, dimensionLessFinalTime, true );

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

            std::string centralBodyOfPropagation = "Sun";

            SelectedAccelerationMap accelerationSettings  = getHaloOrbitAccelerationsMap( );
            basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                        bodyMap, accelerationSettings, { "Spacecraft" }, { centralBodyOfPropagation } );

            std::vector< std::string > centralBodies =  { centralBodyOfPropagation };
            std::vector< std::string > bodiesToPropagate = { "Spacecraft" };

            Eigen::Vector6d initialCartesianState =
                    circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates(
                        secondaryGravitationalParameter, primaryGravitationalParameter,
                        primarySecondaryDistance, currentNormalizedInitialState, dimensionLessInitialTime ) -
                    bodyMap.at( centralBodyOfPropagation )->getEphemeris( )->getCartesianState( initialPropagationTime );

            std::shared_ptr< TranslationalStatePropagatorSettings< double> > propagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate, initialCartesianState,
                      std::make_shared< PropagationTimeTerminationSettings >( finalPropagationTime, true ) );

            std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
                    std::make_shared < numerical_integrators::IntegratorSettings < > >
                    ( numerical_integrators::rungeKutta4, initialPropagationTime, integrationTimeStep );

            SingleArcDynamicsSimulator< > dynamicsSimulator = SingleArcDynamicsSimulator< >(
                        bodyMap, integratorSettings, propagatorSettings );


            std::map< double, Eigen::VectorXd > rawPropagatedStateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

            for( auto stateIterator : rawPropagatedStateHistory )
            {
                propagatedStateHistory[ stateIterator.first ] = stateIterator.second +
                        bodyMap.at( centralBodyOfPropagation )->getEphemeris( )->getCartesianState( stateIterator.first );

                normalizedPropagatedStateHistory[ circular_restricted_three_body_problem::convertDimensionalTimeToDimensionlessTime(
                            stateIterator.first, primaryGravitationalParameter, secondaryGravitationalParameter,
                            primarySecondaryDistance ) ] =
                        circular_restricted_three_body_problem::convertCartesianToCorotatingNormalizedCoordinates(
                            secondaryGravitationalParameter, primaryGravitationalParameter,
                            primarySecondaryDistance, propagatedStateHistory[ stateIterator.first ], stateIterator.first );
            }

            currentNormalizedInitialState = cr3bpStateHistory.rbegin( )->second;

            input_output::writeDataMapToTextFile(
                        unnormalizedCr3bpStateHistory, "cr3bpResultUnnormalized" + std::to_string( i ) +
                        "_" + std::to_string( j ) + ".dat", outputPath );
            input_output::writeDataMapToTextFile(
                        cr3bpStateHistory, "cr3bpResultNormalized.dat" + std::to_string( i ) +
                        "_" + std::to_string( j ) + ".dat", outputPath );
            input_output::writeDataMapToTextFile(
                        propagatedStateHistory, "numericalResultUnnormalized.dat" + std::to_string( i ) +
                        "_" + std::to_string( j ) + ".dat", outputPath );
            input_output::writeDataMapToTextFile(
                        normalizedPropagatedStateHistory, "numericalResultNormalized.dat" + std::to_string( i ) +
                        "_" + std::to_string( j ) + ".dat", outputPath );
        }
    }


    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}


