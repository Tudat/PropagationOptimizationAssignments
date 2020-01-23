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
using namespace tudat;

namespace tudat_applications
{
namespace PropagationOptimization2020
{

enum IntegratorType
{
    Euler,
    RK4,
    RKF45,
    RKF56,
    RKF78,
    RK87,
    ABM,
    BS
};

class ShapeOptimizationProblem
{
public:

    // Constructor
    ShapeOptimizationProblem();

    void setPropagatorIntegrator(TranslationalPropagatorType propagatorType, IntegratorType integratorType);

    std::map< double, Eigen::VectorXd > returnPropagatedStateHistory(std::map< double, Eigen::VectorXd > propagatedStateHistory) const
    {
        return propagatedStateHistory;
    }
    std::map< double, Eigen::VectorXd > returnDependentVariableHistory(std::map< double, Eigen::VectorXd > dependentVariableHistory) const
    {
        return dependentVariableHistory;
    }

    // Fitness function; needs to adhere to Pagmo specifications
    std::vector< double > fitness( std::vector< double >& x ) const;

    // Public member fields
    std::string outputPath;
    // Make mutable so that they can be assigned to in a const class method
    mutable std::map< double, Eigen::VectorXd > propagatedStateHistory;
    mutable std::map< double, Eigen::VectorXd > dependentVariableHistory;

private:

    // Private member fields
    double simulationStartEpoch_;
    Eigen::Vector6d systemInitialState_;
    simulation_setup::NamedBodyMap bodyMap_;
    std::vector< std::string > centralBodies_;
    std::vector< std::string > bodiesToPropagate_;
    basic_astrodynamics::AccelerationMap accelerationModelMap_;
    std::shared_ptr< IntegratorSettings< > > integratorSettings_;
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings_;

};

}
}
