#include "Marmot/GradientVonMises.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include <Eigen/Dense>
#include <functional>
#include <vector>

using namespace Marmot;
using namespace Marmot::Materials;
using namespace Marmot::Testing;

using Mat = MarmotMaterialGradientPlasticityHypoElastic< 1 >;
using Res = Mat::response;
using Tan = Mat::tangents;
using Inc = Mat::increment;

static GradientVonMises makeMaterial( const std::vector< double >& props )
{
  return GradientVonMises( props.data(), static_cast< int >( props.size() ), 1 );
}

static void initializeResponse( const GradientVonMises& mat, Res& res, std::vector< double >& stateVars )
{
  stateVars.assign( mat.getNumberOfRequiredStateVars(), 0.0 );
  res.stress               = Vector6d::Zero();
  res.f                    = Eigen::Vector< double, 1 >::Zero();
  res.stateVars            = stateVars.data();
  res.elasticEnergyDensity = 0.0;
  res.dissipation          = 0.0;
}

static std::tuple< Vector6d, double, Tan > evaluateModel( const GradientVonMises& mat,
                                                          const Vector6d&         dStrain,
                                                          const double            dLambda )
{
  std::vector< double > stateVars;
  Res                   res;
  Tan                   tan;
  Inc                   inc;

  initializeResponse( mat, res, stateVars );

  inc.dStrain             = dStrain;
  inc.dLambda( 0 )        = dLambda;
  inc.laplaceDLambda( 0 ) = 0.0;
  inc.time                = 0.0;
  inc.dT                  = 1.0;

  mat.computeStress( res, tan, inc );

  return { res.stress, res.f( 0 ), tan };
}

void testGradientVonMisesTangentByNumericalDifferentiation()
{
  const std::vector< double > props = { 210000.0, 0.3, 250.0, 1500.0, 25.0, 7850.0, 0.0 };
  const auto                  mat   = makeMaterial( props );

  Vector6d dStrain     = Vector6d::Zero();
  dStrain( 0 )         = 1.0e-6;
  dStrain( 1 )         = -5.0e-7;
  dStrain( 2 )         = -5.0e-7;
  dStrain( 3 )         = 2.0e-7;
  const double dLambda = 0.0;

  auto [stressRef, fRef, tanRef] = evaluateModel( mat, dStrain, dLambda );

  using namespace Marmot::NumericalAlgorithms::Differentiation;
  const vector_to_vector_function_type residualFunction = [&]( const Eigen::VectorXd& x ) {
    Vector6d     strainEval = x.head< 6 >();
    const double lambdaEval = x( 6 );
    auto [stress, f, tan]   = evaluateModel( mat, strainEval, lambdaEval );

    Eigen::VectorXd y( 7 );
    y.head< 6 >() = stress;
    y( 6 )        = f;
    return y;
  };

  Eigen::VectorXd x( 7 );
  x.head< 6 >() = dStrain;
  x( 6 )        = dLambda;

  const Eigen::MatrixXd jacNum = centralDifference( residualFunction, x );

  const Matrix6d dStressddStrainNum = jacNum.block< 6, 6 >( 0, 0 );
  const Vector6d dStressddLambdaNum = jacNum.block< 6, 1 >( 0, 6 );

  throwExceptionOnFailure( checkIfEqual< double >( tanRef.dStressddStrain, dStressddStrainNum, 5e-5 ),
                           "Mismatch in dStressddStrain vs numerical Jacobian in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  throwExceptionOnFailure( checkIfEqual< double >( tanRef.dStressddLambda.col( 0 ), dStressddLambdaNum, 5e-5 ),
                           "Mismatch in dStressddLambda vs numerical Jacobian in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  throwExceptionOnFailure( std::isfinite( stressRef.norm() ) && std::isfinite( fRef ),
                           "Non-finite reference response encountered in " + std::string( __PRETTY_FUNCTION__ ) );
}

void testGradientVonMisesPlasticTangentByNumericalDifferentiation()
{
  const std::vector< double > props = { 210000.0, 0.3, 250.0, 1500.0, 25.0, 7850.0, 0.0 };
  const auto                  mat   = makeMaterial( props );

  Vector6d dStrain     = Vector6d::Zero();
  dStrain( 0 )         = 3.0e-3;
  dStrain( 1 )         = -1.5e-3;
  dStrain( 2 )         = -1.5e-3;
  dStrain( 3 )         = 1.2e-3;
  const double dLambda = 2.0e-4;

  auto [stressRef, fRef, tanRef] = evaluateModel( mat, dStrain, dLambda );

  throwExceptionOnFailure( fRef > 0.0,
                           "Chosen plastic test point is not yielding in " + std::string( __PRETTY_FUNCTION__ ) );

  using namespace Marmot::NumericalAlgorithms::Differentiation;
  const vector_to_vector_function_type stressOfStrain = [&]( const Eigen::VectorXd& x ) {
    const Vector6d strainEval = x.head< 6 >();
    auto [stress, f, tan]     = evaluateModel( mat, strainEval, dLambda );
    return stress;
  };

  Eigen::VectorXd x( 6 );
  x                                  = dStrain;
  const Eigen::MatrixXd jacStrainNum = centralDifference( stressOfStrain, x );

  const vector_to_vector_function_type stressOfLambda = [&]( const Eigen::VectorXd& xLambda ) {
    auto [stress, f, tan] = evaluateModel( mat, dStrain, xLambda( 0 ) );
    return stress;
  };

  Eigen::VectorXd xLambda( 1 );
  xLambda( 0 )                       = dLambda;
  const Eigen::MatrixXd jacLambdaNum = centralDifference( stressOfLambda, xLambda );

  const double relErrStrain = ( tanRef.dStressddStrain - jacStrainNum ).norm() /
                              std::max( tanRef.dStressddStrain.norm(), 1.0 );
  const double relErrLambda = ( tanRef.dStressddLambda.col( 0 ) - jacLambdaNum.col( 0 ) ).norm() /
                              std::max( tanRef.dStressddLambda.col( 0 ).norm(), 1.0 );

  throwExceptionOnFailure( relErrStrain < 5e-2,
                           MakeString() << "Mismatch in plastic dStressddStrain vs numerical Jacobian (rel err = "
                                        << relErrStrain << ") in " << __PRETTY_FUNCTION__ );

  throwExceptionOnFailure( relErrLambda < 5e-2,
                           MakeString() << "Mismatch in plastic dStressddLambda vs numerical Jacobian (rel err = "
                                        << relErrLambda << ") in " << __PRETTY_FUNCTION__ );

  throwExceptionOnFailure( std::isfinite( stressRef.norm() ),
                           "Non-finite plastic-point stress encountered in " + std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  const std::vector< std::function< void() > > tests = {
    testGradientVonMisesTangentByNumericalDifferentiation,
    testGradientVonMisesPlasticTangentByNumericalDifferentiation,
  };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
