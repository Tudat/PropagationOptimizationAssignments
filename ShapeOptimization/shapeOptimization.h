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
//std::shared_ptr< HypersonicLocalInclinationAnalysis > getCapsuleCoefficientInterface(
//        const std::shared_ptr< geometric_shapes::Capsule > capsule,
//        const std::string directory,
//        const std::string filePrefix,
//        const bool useNewtonianMethodForAllPanels = true );

void setVehicleShapeParameters(
        std::vector< double > shapeParameters,
        const NamedBodyMap& bodyMap );

//! Class to set the aerodynamic angles of the capsule (default: all angles 0)
class CapsuleAerodynamicGuidance: public aerodynamics::AerodynamicGuidance
{
public:

    //! Constructor
    CapsuleAerodynamicGuidance(
            const NamedBodyMap bodyMap_,
            const double fixedAngleOfAttack ):bodyMap_( bodyMap_ ), fixedAngleOfAttack_( fixedAngleOfAttack )
    {

    }

    //! The aerodynamic angles are to be computed here
    void updateGuidance( const double time )
    {
        currentAngleOfAttack_ = fixedAngleOfAttack_;
        currentAngleOfSideslip_ = 0.0;
        currentBankAngle_ = 0.0;

    }

private:

    //! List of body objects that constitute the environment
    NamedBodyMap bodyMap_;

    //! Fixed angle of attack that is to be used by vehicle
    double fixedAngleOfAttack_;
};



class ShapeOptimizationProblem
{
public:

    // Constructor
    ShapeOptimizationProblem(
            const simulation_setup::NamedBodyMap bodyMap,
            const std::shared_ptr< IntegratorSettings< > > integratorSettings,
            const std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings );

    ShapeOptimizationProblem( );

    std::map< double, Eigen::VectorXd > getLastRunPropagatedStateHistory( ) const
    {
        return propagatedStateHistory;
    }
    std::map< double, Eigen::VectorXd > getLastRunDependentVariableHistory( ) const
    {
        return dependentVariableHistory;
    }

    // Fitness function; needs to adhere to Pagmo specifications
    std::vector< double > fitness( std::vector< double >& x ) const;



private:

    simulation_setup::NamedBodyMap bodyMap_;
    std::shared_ptr< IntegratorSettings< > > integratorSettings_;
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings_;

    // Make mutable so that they can be assigned to in a const class method
    mutable std::map< double, Eigen::VectorXd > propagatedStateHistory;
    mutable std::map< double, Eigen::VectorXd > dependentVariableHistory;

};

}
}
