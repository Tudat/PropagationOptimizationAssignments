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

int main()
{
    // Construct problem
    ShapeOptimizationProblem prob{};

    // Loop for the integrator/propagator settings
    const unsigned int numberPropagatorCases = 7;
    const unsigned int numberIntegratorCases = 3;

    TranslationalPropagatorType propagatorType;
    IntegratorType integratorType;

    // Perform benchmark
    prob.setPropagatorIntegrator(cowell, RKF78);
    input_output::writeDataMapToTextFile( prob.propagatedStateHistory, "stateHistory_benchmark.dat", prob.outputPath );
    input_output::writeDataMapToTextFile( prob.dependentVariableHistory, "dependentVariables_benchmark.dat", prob.outputPath );

    // Decision vector is required by Pagmo specification. Since we are calling the fitness function directly,
    // pass it a dummy decision vector {0.0}
    std::vector< double > decisionVector{0.0};
    prob.fitness( decisionVector );

    // Run different combinations of propagator/integrator settings
    for( unsigned int i = 0; i < numberPropagatorCases; i++ )
    {
        for( unsigned int j = 0; j < numberIntegratorCases; j++ )
        {

            switch( i )
            {
                case 0:
                    propagatorType = cowell;
                    break;
                case 1:
                    propagatorType = encke;
                    break;
                case 2:
                    propagatorType = gauss_keplerian;
                    break;
                case 3:
                    propagatorType = gauss_modified_equinoctial;
                    break;
                case 4:
                    propagatorType = unified_state_model_exponential_map;
                    break;
                case 5:
                    propagatorType = unified_state_model_modified_rodrigues_parameters;
                    break;
                case 6:
                    propagatorType = unified_state_model_quaternions;
                    break;
                default:
                    propagatorType = cowell;
                    break;
            }

            switch( j )
            {
                case 0:
                    integratorType = RK4;
                    break;
                case 1:
                    integratorType = RKF78;
                    break;
                case 2:
                    integratorType = ABM;
                    break;
                default:
                    integratorType = RK4;
                    break;
            }

            // Set integrator & propagator settings
            prob.setPropagatorIntegrator(propagatorType, integratorType);

            // Propagate trajectory using defined settings
            prob.fitness( decisionVector );

            // TODO: retrieve state & dependent variable history from fitness function (const function!)
            input_output::writeDataMapToTextFile( prob.propagatedStateHistory, "stateHistory_" + std::to_string(i) + "_" + std::to_string(j) + ".dat", prob.outputPath );
            input_output::writeDataMapToTextFile( prob.dependentVariableHistory, "dependentVariables_" + std::to_string(i) + ".dat", prob.outputPath );
        }
    }
}
