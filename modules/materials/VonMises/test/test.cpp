#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/VonMises.h"

using namespace Marmot::Testing;
using namespace Marmot::Solvers;

void testVonMisesCoordinateInvariance()
{
  // material properties
  std::vector< double > materialProperties = { 210000., 0.3, 200., 2100., 20., 20 };
  auto                  solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string           matName            = "VONMISES";

  // create material point solver instance
  auto solver = MarmotMaterialPointSolverHypoElastic( matName,
                                                      &materialProperties[0],
                                                      materialProperties.size(),
                                                      solveropts );

  // define a step
  MarmotMaterialPointSolverHypoElastic::Step step;

  // define step parameters
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.01, 0., 0., 0., 0.03, 0.0 };
  step.dTStart                     = 1.0;

  // add step to solver
  solver.addStep( step );
  // check coordinate invariance with Turbokreisel
  throwExceptionOnFailure( spinTurbokreisel( solver, 1e-10, 1e-8 ), "Turbokreisel failed!" );
}

void testVonMises()
{
  // material properties
  std::vector< double > materialProperties = { 210000., 0.3, 200., 2100., 20., 20 };
  auto                  solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string           matName            = "VONMISES";

  // create material point solver instance
  auto solver = MarmotMaterialPointSolverHypoElastic( matName,
                                                      &materialProperties[0],
                                                      materialProperties.size(),
                                                      solveropts );

  // define a step
  MarmotMaterialPointSolverHypoElastic::Step step;

  // define step parameters
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.00839244, 0.00089344, -0.00703916, 0.00013635, 0.00160548, 0.00572825 };
  step.dTStart                     = 1.0;

  // add step to solver
  solver.addStep( step );

  solver.solve();

  // read history
  auto history = solver.getHistory();
  // reference solution
  const Marmot::Vector6d stressTarget =
    { 511.262747695, 395.408929527, 272.856322779, 1.0532516407, 12.4017194288, 44.2485420671 };

  // compare solutions
  throwExceptionOnFailure( checkIfEqual< double >( history.back().stress, stressTarget, 1e-9 ),
                           "comparison with reference solution failed" );
}

void testVonMisesDensityAndDamping()
{
  const std::vector< double > materialProperties = { 210000., 0.3, 200., 2100., 20., 20., 7.5, 0.15 };
  auto                        mat = Marmot::Materials::VonMisesModel( materialProperties.data(),
                                                materialProperties.size(),
                                                1 );

  throwExceptionOnFailure( checkIfEqual( mat.getDensity(), 7.5, 1e-12 ),
                           "density retrieval failed for VonMises in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual( mat.getDampingCoefficient(), 0.15, 1e-12 ),
                           "damping retrieval failed for VonMises in " + std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  std::vector< std::function< void( void ) > > tests = {
    testVonMises,
    testVonMisesCoordinateInvariance,
    testVonMisesDensityAndDamping,
  };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
