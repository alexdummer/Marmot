#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "Marmot/MarmotMaterialPointSolverFiniteStrain.h"
#include "fastor_nanobind_caster.h"

namespace nb = nanobind;
using Solver = Marmot::Solvers::MarmotMaterialPointSolverFiniteStrain;

void bind_finite_strain_solver( nb::module_& m )
{
  nb::class_< Solver > solver( m, "FiniteStrainSolver" );

  nb::class_< Solver::Step >( solver, "Step" )
    .def( nb::init<>() )
    .def_rw( "gradUIncrementTarget", &Solver::Step::gradUIncrementTarget )
    .def_rw( "stressIncrementTarget", &Solver::Step::stressIncrementTarget )
    .def_rw( "isGradUComponentControlled", &Solver::Step::isGradUComponentControlled )
    .def_rw( "isStressComponentControlled", &Solver::Step::isStressComponentControlled )
    .def_rw( "timeStart", &Solver::Step::timeStart )
    .def_rw( "timeEnd", &Solver::Step::timeEnd )
    .def_rw( "dTStart", &Solver::Step::dTStart )
    .def_rw( "dTMin", &Solver::Step::dTMin )
    .def_rw( "dTMax", &Solver::Step::dTMax )
    .def_rw( "maxIncrements", &Solver::Step::maxIncrements )
    .def( "checkControl", &Solver::Step::checkControl );

  nb::class_< Solver::Increment >( solver, "Increment" )
    .def( nb::init<>() )
    .def_rw( "gradUIncrement", &Solver::Increment::gradUIncrement )
    .def_rw( "stressIncrement", &Solver::Increment::stressIncrement )
    .def_rw( "isGradUComponentControlled", &Solver::Increment::isGradUComponentControlled )
    .def_rw( "isStressComponentControlled", &Solver::Increment::isStressComponentControlled )
    .def_rw( "timeOld", &Solver::Increment::timeOld )
    .def_rw( "dT", &Solver::Increment::dT );

  nb::class_< Solver::SolverOptions >( solver, "SolverOptions" )
    .def( nb::init<>() )
    .def_rw( "maxIterations", &Solver::SolverOptions::maxIterations )
    .def_rw( "residualTolerance", &Solver::SolverOptions::residualTolerance )
    .def_rw( "correctionTolerance", &Solver::SolverOptions::correctionTolerance );

  nb::class_< Solver::HistoryEntry >( solver, "HistoryEntry" )
    .def( nb::init<>() )
    .def_ro( "time", &Solver::HistoryEntry::time )
    .def_ro( "stress", &Solver::HistoryEntry::stress )
    .def_ro( "F", &Solver::HistoryEntry::F )
    .def_ro( "dTau_dF", &Solver::HistoryEntry::dTau_dF )
    .def_ro( "stateVars", &Solver::HistoryEntry::stateVars )
    .def( "print", &Solver::HistoryEntry::print );

  solver
    .def( "__init__",
          []( Solver*                                            t,
              const std::string&                                 name,
              nb::ndarray< double, nb::ndim< 1 >, nb::c_contig > props,
              const Solver::SolverOptions&                       opts ) {
            new ( t ) Solver( name, props.data(), static_cast< int >( props.size() ), opts );
          } )
    .def( "addStep", &Solver::addStep )
    .def( "getSteps", &Solver::getSteps )
    .def( "clearSteps", &Solver::clearSteps )
    .def( "setInitialState", &Solver::setInitialState )
    .def( "getNumberOfStateVariables", &Solver::getNumberOfStateVariables )
    .def( "resetToInitialState", &Solver::resetToInitialState )
    .def( "solve", &Solver::solve )
    .def( "getHistory", &Solver::getHistory )
    .def( "clearHistory", &Solver::clearHistory )
    .def( "printHistory", &Solver::printHistory )
    .def( "exportHistoryToCSV", &Solver::exportHistoryToCSV );
}
