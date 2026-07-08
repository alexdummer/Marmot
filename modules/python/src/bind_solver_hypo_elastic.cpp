#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"

namespace nb = nanobind;
using Solver = Marmot::Solvers::MarmotMaterialPointSolverHypoElastic;

void bind_hypo_elastic_solver( nb::module_& m )
{
  nb::class_< Solver > solver( m, "HypoElasticSolver" );

  nb::class_< Solver::Step >( solver, "Step" )
    .def( nb::init<>() )
    .def_rw( "strainIncrementTarget", &Solver::Step::strainIncrementTarget )
    .def_rw( "stressIncrementTarget", &Solver::Step::stressIncrementTarget )
    .def_rw( "isStrainComponentControlled", &Solver::Step::isStrainComponentControlled )
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
    .def_rw( "strainIncrement", &Solver::Increment::strainIncrement )
    .def_rw( "stressIncrement", &Solver::Increment::stressIncrement )
    .def_rw( "isStrainComponentControlled", &Solver::Increment::isStrainComponentControlled )
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
    .def_ro( "strain", &Solver::HistoryEntry::strain )
    .def_ro( "dStressdStrain", &Solver::HistoryEntry::dStressdStrain )
    .def_ro( "stateVars", &Solver::HistoryEntry::stateVars )
    .def( "print", &Solver::HistoryEntry::print );

  solver
    .def( "__init__",
          []( Solver*                             t,
              std::string                         name,
              nb::ndarray< double, nb::c_contig > props,
              const Solver::SolverOptions& opts ) { new ( t ) Solver( name, props.data(), props.size(), opts ); } )
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
