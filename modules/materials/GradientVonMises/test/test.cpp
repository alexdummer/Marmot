#include "Marmot/GradientVonMises.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include <vector>

using namespace Marmot::Testing;
using namespace Marmot::Materials;
using namespace Marmot;

using Mat = MarmotMaterialGradientPlasticityHypoElastic;
using Res = Mat::response;
using Tan = Mat::tangents;
using Inc = Mat::increment;

void testElasticStep()
{
  const std::vector< double > props = { 210000., 0.3, 400., -2000., 0.05, 7850.0, 0.01 };
  GradientVonMises            mat( props.data(), props.size(), 1 );
  std::vector< double >       stateVars( mat.getNumberOfRequiredStateVars(), 0.0 );

  Res res;
  Tan tan;
  Inc inc;
  inc.dStrain      = Vector6d::Zero();
  inc.dStrain( 0 ) = 1e-4;
  inc.K( 0 )       = 0.1;
  inc.dK( 0 )      = 0.0;
  inc.time         = 0.0;
  inc.dT           = 1.0;

  res.stress    = Vector6d::Zero();
  res.KLocal    = Eigen::Vector< double, 1 >::Zero();
  res.c         = Eigen::Vector< double, 1 >::Zero();
  res.stateVars = stateVars.data();

  mat.computeStress( res, tan, inc );

  const Matrix6d C = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( props[0], props[1] );
  throwExceptionOnFailure( checkIfEqual< double >( res.stress, C * inc.dStrain, 1e-10 ), "elastic stress mismatch" );
  throwExceptionOnFailure( checkIfEqual( res.KLocal( 0 ), 0.0, 1e-14 ), "KLocal should remain zero in elastic step" );
  throwExceptionOnFailure( checkIfEqual( res.c( 0 ), props[4] * props[4], 1e-14 ), "wrong gradient parameter c" );
}

void testPlasticStepAndCoupling()
{
  const std::vector< double > props = { 210000., 0.3, 200., -2500., 0.05, 7850.0, 0.01 };
  GradientVonMises            mat( props.data(), props.size(), 1 );
  std::vector< double >       stateVars( mat.getNumberOfRequiredStateVars(), 0.0 );

  Res res;
  Tan tan;
  Inc inc;
  inc.dStrain      = Vector6d::Zero();
  inc.dStrain( 0 ) = 2e-3;
  inc.K( 0 )       = 0.2;
  inc.dK( 0 )      = 0.0;
  inc.time         = 0.0;
  inc.dT           = 1.0;

  res.stress    = Vector6d::Zero();
  res.KLocal    = Eigen::Vector< double, 1 >::Zero();
  res.c         = Eigen::Vector< double, 1 >::Zero();
  res.stateVars = stateVars.data();

  mat.computeStress( res, tan, inc );

  throwExceptionOnFailure( res.KLocal( 0 ) > 0.0, "KLocal should increase in plastic step" );
  throwExceptionOnFailure( tan.dStressddK.col( 0 ).norm() > 0.0, "stress/nonlocal coupling should be active" );
  throwExceptionOnFailure( tan.dKLocalddStrain.norm() > 0.0, "KLocal/strain coupling should be active" );
}

int main()
{
  std::vector< std::function< void() > > tests = { testElasticStep, testPlasticStepAndCoupling };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
