#include "Marmot/WCPPlasticity.h"
#include "Marmot/MarmotAutomaticDifferentiation.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/NewtonConvergenceChecker.h"

namespace Marmot {

  struct FailedToConverge : std::exception {};
} // namespace Marmot

namespace Marmot::Materials {

  using namespace Eigen;
  WCPPlasticity::WCPPlasticity( const MaterialParameters& matParams ) : mp( matParams ) {}

  WCPPlasticity::ReturnMapResult WCPPlasticity::performReturnMapping( const MaterialState& trialState,
                                                                      const MaterialState& initialGuessState,
                                                                      const Matrix6d&      elasticStiffness,
                                                                      const Matrix6d&      elasticCompliance )
  {

    YieldSurfFlagArr activeSurfaces = evaluateYieldFunctions( trialState );

    return returnMappingAttempt( trialState, initialGuessState, elasticStiffness, elasticCompliance, activeSurfaces );
  }

  WCPPlasticity::ReturnMapResult WCPPlasticity::performSmartReturnMapping( const MaterialState& trialState,
                                                                           const MaterialState& oldState,
                                                                           const Matrix6d&      elasticStiffness,
                                                                           const Matrix6d&      elasticCompliance )
  {
    try {
      return performReturnMapping( trialState, trialState, elasticStiffness, elasticCompliance );
    }
    catch ( const ReturnMappingFailedException& e ) {
      return performReturnMapping( trialState, oldState, elasticStiffness, elasticCompliance );
    }
  }

  WCPPlasticity::ReturnMapResult WCPPlasticity::retryReturnMapping( const MaterialState& trialState,
                                                                    const MaterialState& initialGuessState,
                                                                    const Matrix6d&      elasticStiffness,
                                                                    const Matrix6d&      elasticCompliance )
  {

    /* YieldSurfFlagArr newActiveSurfaces; */
    /* if ( !yieldSurfCombiManager.getAnotherYieldFlagCombination( newActiveSurfaces ) ) */
    /*   throw ReturnMappingFailedException(); */

    /* return returnMappingAttempt( trialState, */
    /*                              initialGuessState, */
    /*                              elasticStiffness, */
    /*                              elasticCompliance, */
    /*                              newActiveSurfaces ); */
    return ReturnMapResult();
  }

  WCPPlasticity::ReturnMapResult WCPPlasticity::returnMappingAttempt( const MaterialState& trialState,
                                                                      const MaterialState& initialGuessState,
                                                                      const Matrix6d&      elasticStiffness,
                                                                      const Matrix6d&      elasticCompliance,
                                                                      YieldSurfFlagArr&    activeSurfaces )
  {
    /* std::cout << "activeSurfaces" << activeSurfaces << std::endl; */
    const auto nActiveYieldSurfaces = activeSurfaces.count();
    const auto sizeEq               = nStress + nActiveYieldSurfaces * 2;

    /* yieldSurfCombiManager.markYieldFlagCombinationAsUsed( activeSurfaces ); */

    // solution vector X, intialize with initial guess
    Eigen::VectorXd X = VectorXd::Zero( sizeEq );
    X.head( nStress ) = initialGuessState.stress;

    Eigen::VectorXd dX = VectorXd::Zero( sizeEq );
    // aux vectors for solving the nonlinear equation system using Newton
    Eigen::VectorXd R( sizeEq );
    Eigen::MatrixXd dR_dX_( sizeEq, sizeEq );

    std::tie( R, dR_dX_ ) = dR_dX( X, trialState, elasticStiffness, activeSurfaces );

    NumericalAlgorithms::NewtonConvergenceChecker newtonChecker( createResidualScaleVector( R ),
                                                                 10,    // nMaxInnerNewtonCycles,
                                                                 15,    // nMaxInnerNewtonCyclesAlt,
                                                                 1e-12, // innerNewtonTol,
                                                                 1e-10, // innerNewtonRTol,
                                                                 1e-8,  // innerNewtonTolAlt,
                                                                 1e-6   // innerNewtonRTolAlt
    );

    /* std::cout << "R: " << std::endl << R.transpose() << std::endl; */
    /* std::cout << "dR_dX_: "<< std::endl << dR_dX_ << std::endl; */

    /* std::cout << "R: " << std::endl << R.transpose() << std::endl; */
    /* std::cout << "dR_dX_: "<< std::endl << dR_dX_ << std::endl; */

    try {

      int iterationCounter = 0;
      while ( newtonChecker.isConverged( R, X, dX, iterationCounter ) == false ) {

        if ( iterationCounter > 15 )
          throw FailedToConverge();

        if ( Math::isNaN( R ) )
          throw FailedToConverge();

        if ( Math::isNaN( dR_dX_ ) )
          throw FailedToConverge();

        dX = -dR_dX_.colPivHouseholderQr().solve( R );
        X += dX;
        std::tie( R, dR_dX_ ) = dR_dX( X, trialState, elasticStiffness, activeSurfaces );
        iterationCounter++;
      }
    }

    catch ( const FailedToConverge& e ) {
      throw ReturnMappingFailedException();
      /* std::exit(1); */
      /* return retryReturnMapping( trialState, initialGuessState, elasticStiffness, elasticCompliance ); */
    }

    /* if ( ( X.tail( nActiveYieldSurfaces ).array() < 1e-5 ).any() ) { */
    /*   return retryReturnMapping( trialState, initialGuessState, elasticStiffness, elasticCompliance ); */
    /* } */

    const MaterialState newMaterialState = extractMaterialState( X, trialState, activeSurfaces );

    /* std::cout << "newMaterialState.stress: " << std::endl << newMaterialState.stress.transpose() << std::endl; */

    /* std::exit(1); */

    /* if ( checkIfYielding( newMaterialState ) ) { */
    /*   return retryReturnMapping( trialState, initialGuessState, elasticStiffness, elasticCompliance ); */
    /* } */

    return ReturnMapResult{ .materialState     = newMaterialState,
                            .dStrainPlastic    = computePlasticStrainIncrement( newMaterialState.stress,
                                                                             trialState.stress,
                                                                             elasticCompliance ),
                            .algorithmicModuli = computeAlgorithmicModuli( dR_dX_,
                                                                           elasticStiffness,
                                                                           elasticCompliance ) };
  } // namespace Marmot::Materials
  VectorXd WCPPlasticity::createResidualScaleVector( const VectorXd& LHS )
  {
    const auto nRows = LHS.rows();
    VectorXd   S( nRows );
    S.setOnes();
    /* const double maxAbsStress = LHS.segment< nStress >( idxS ).array().abs().maxCoeff(); */
    /* S.segment< nStress >( idxS ).setConstant( 1. / std::max( 1.0, maxAbsStress ) ); */
    /* S( idxInternal ) = 1. / std::max( 1.0, std::abs( LHS( idxInternal ) ) ); */

    /* const int nActiveYieldSurfaces = nRows - nStress - nInternalVariables; */

    // yielding is checked anyway to a certain toS.setConstant( 1. );

    return S;
  }

  WCPPlasticity::AlgorithmicModuli WCPPlasticity::computeAlgorithmicModuli( const MatrixXd& dF_dX,
                                                                            const Matrix6d& elasticStiffness,
                                                                            const Matrix6d& elasticCompliance )
  {

    const auto&       Cel = elasticStiffness;
    AlgorithmicModuli result;

    const auto      sizeEq = dF_dX.rows();
    Eigen::MatrixXd dYdE   = Eigen::MatrixXd::Identity( sizeEq, sizeEq );

    const Eigen::MatrixXd dXdE = dF_dX.colPivHouseholderQr().solve( dYdE );
    /* std::cout << "dFdX: " << std::endl << dF_dX << std::endl; */
    /* std::cout << "dXdE: " << std::endl << dXdE.transpose() << std::endl; */
    result.dStress_dStrain        = dXdE.block< nStress, nStress >( idxS, idxS ) * Cel;
    result.dStrainPlastic_dStrain = -elasticCompliance * result.dStress_dStrain + Matrix6d::Identity();
    /* std::cout << "dStress_dStrain: " << std::endl << result.dStress_dStrain << std::endl; */
    /* std::exit(1); */
    return result;
  }
  Vector6d WCPPlasticity::computePlasticStrainIncrement( const Vector6d& stressNew,
                                                         const Vector6d& trialStress,
                                                         const Matrix6d& elasticCompliance )
  {
    return elasticCompliance * ( trialStress - stressNew );
  }

  WCPPlasticity::MaterialState WCPPlasticity::extractMaterialState( const Eigen::VectorXd&  X,
                                                                    const MaterialState&    trialState,
                                                                    const YieldSurfFlagArr& activeSurfaces )
  {
    MaterialState newMaterialState = trialState;
    /* std::cout << "X: " << std::endl << X.transpose() << std::endl; */
    newMaterialState.stress = X.segment< nStress >( idxS );

    int idxInternal_ = idxInternal;
    if ( activeSurfaces( 0 ) ) {
      newMaterialState.alphaR = X( idxInternal_ );
      idxInternal_ += 2;
    }
    if ( activeSurfaces( 1 ) ) {
      newMaterialState.alphaT = X( idxInternal_ );
      idxInternal_ += 2;
    }
    if ( activeSurfaces( 2 ) ) {
      newMaterialState.alphaT = X( idxInternal_ );
      idxInternal_ += 2;
    }
    return newMaterialState;
  };
  std::tuple< VectorXd, MatrixXd > WCPPlasticity::dR_dX( const VectorXd&                  X,
                                                         const MaterialState&             trialState,
                                                         const Matrix6d&                  elasticStiffness,
                                                         WCPPlasticity::YieldSurfFlagArr& activeSurfaces )
  {

    using namespace Marmot::AutomaticDifferentiation;

    vector_to_vector_function_type_dual func = [&]( const autodiff::VectorXdual& X_ ) {
      return R( X_, trialState, elasticStiffness, activeSurfaces );
    };

    return AutomaticDifferentiation::jacobian( func, X );
  }
} // namespace Marmot::Materials
