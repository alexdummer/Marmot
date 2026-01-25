#include "Marmot/UnistrandPlasticity.h"
#include "Marmot/MarmotAutomaticDifferentiation.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/NewtonConvergenceChecker.h"

namespace Marmot {

  struct FailedToConverge : std::exception {};
} // namespace Marmot

namespace Marmot::Materials {

  using namespace Eigen;
  UnistrandPlasticity::UnistrandPlasticity( const strength& st_, const hardening& hard, const double l ) : st( st_ )
  {

    // get compressive exp parameters direct from hardening law
    mu_.mu1c = hard.mu1c;
    mu_.mu2c = hard.mu2c;
    mu_.mu3c = hard.mu3c;

    // convert tensile fracture energies to exp parameters
    mu_.mu1t = -l * st.r11t / hard.Gf1 * 0;
    mu_.mu2t = -l * st.r22t / hard.Gf2 * 0;
    mu_.mu3t = -l * st.r33t / hard.Gf3 * 0;

    // convert shear fracture energies to exp parameters
    mu_.mu12 = -l * st.s12 / hard.Gs12 * 0;
    mu_.mu13 = -l * st.s13 / hard.Gs13 * 0;
    mu_.mu23 = -l * st.s23 / hard.Gs23 * 0;
  }

  UnistrandPlasticity::ReturnMapResult UnistrandPlasticity::performReturnMapping(
    const MaterialState& trialState,
    const MaterialState& initialGuessState,
    const Matrix6d&      elasticStiffness,
    const Matrix6d&      elasticCompliance )
  {

    YieldSurfFlagArr activeSurfaces = evaluateYieldFunctions( trialState );

    return returnMappingAttempt( trialState, initialGuessState, elasticStiffness, elasticCompliance, activeSurfaces );
  }

  UnistrandPlasticity::ReturnMapResult UnistrandPlasticity::performSmartReturnMapping(
    const MaterialState& trialState,
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

  UnistrandPlasticity::ReturnMapResult UnistrandPlasticity::returnMappingAttempt(
    const MaterialState& trialState,
    const MaterialState& initialGuessState,
    const Matrix6d&      elasticStiffness,
    const Matrix6d&      elasticCompliance,
    YieldSurfFlagArr&    activeSurfaces )
  {

    // const auto nActiveYieldSurfaces = activeSurfaces.count();
    const auto sizeEq               = nStress + 9 * 2;

    yieldSurfCombiManager.markYieldFlagCombinationAsUsed( activeSurfaces );

    // solution vector X, intialize with initial guess
    Eigen::VectorXd X = VectorXd::Zero( sizeEq );
    X.head( nStress ) = initialGuessState.stress;

    int idxCurrentInternal = idxInternal;
    for ( int i = 0; i < nYieldSurfaces; i++ ) {
      // if ( activeSurfaces( i ) ) {
      X( idxCurrentInternal ) = initialGuessState.alpha( i );
      idxCurrentInternal += 2;
      // } 
    }

    Eigen::VectorXd dX = VectorXd::Zero( sizeEq );
    // aux vectors for solving the nonlinear equation system using Newton
    Eigen::VectorXd R( sizeEq );
    Eigen::MatrixXd dR_dX_( sizeEq, sizeEq );

    std::tie( R, dR_dX_ ) = dR_dX( X, trialState, elasticStiffness, activeSurfaces );

    NumericalAlgorithms::NewtonConvergenceChecker newtonChecker( createResidualScaleVector( R ),
                                                                 nMaxInnerNewtonCycles,
                                                                 nMaxInnerNewtonCyclesAlt,
                                                                 innerNewtonTol,
                                                                 innerNewtonRTol,
                                                                 innerNewtonTolAlt,
                                                                 innerNewtonRTolAlt );

    try {

      int iterationCounter = 0;
      while ( newtonChecker.isConverged( R, X, dX, iterationCounter ) == false ) {

        if ( iterationCounter > nMaxInnerNewtonCyclesAlt ){
          // std::cout << "UnistrandPlasticity::returnMappingAttempt: Failed to converge within "
          //              << nMaxInnerNewtonCyclesAlt << " iterations." <<std::endl;
          //   std::cout << "R: " << R.transpose() <<std::endl;
          //   std::cout << "X: " << X.transpose() <<std::endl;
          //   std::cout << "dX: " << dX.transpose() <<std::endl;
          //   std::cout << "Iteration: " << iterationCounter <<std::endl;
          //   std::cout << "dR_dX_: " << dR_dX_ <<std::endl;
            throw FailedToConverge();
        }

        if ( Math::isNaN( R ) ){
            // std::cout << "UnistrandPlasticity::returnMappingAttempt: Residual contains NaN." <<std::endl;
            // std::cout << "R: " << R.transpose() <<std::endl;
            // std::cout << "X: " << X.transpose() <<std::endl;
            // std::cout << "dX: " << dX.transpose() <<std::endl;
            // std::cout << "Iteration: " << iterationCounter <<std::endl;
            // std::cout << "dR_dX_: " << dR_dX_ <<std::endl;
            throw FailedToConverge();
        }
        if ( Math::isNaN( dR_dX_ ) ){
            // std::cout << "UnistrandPlasticity::returnMappingAttempt: Jacobian contains NaN." <<std::endl;
            // std::cout << "R: " << R.transpose() <<std::endl;
            // std::cout << "X: " << X.transpose() <<std::endl;
            // std::cout << "dX: " << dX.transpose() <<std::endl;
            // std::cout << "Iteration: " << iterationCounter <<std::endl;
            // std::cout << "dR_dX_: " << dR_dX_ <<std::endl;
          throw FailedToConverge();}

        dX = -dR_dX_.colPivHouseholderQr().solve( R );

        for ( int i = 0; i < nYieldSurfaces; i++ ) {
            // if the yield surface is not active, do not update the corresponding internal variables
            const int idxAlpha = idxInternal + i * 2;
            dX( idxAlpha )     = std::max( 0.0, dX( idxAlpha ) );
          }
        
        X += dX;
        std::tie( R, dR_dX_ ) = dR_dX( X, trialState, elasticStiffness, activeSurfaces );
        iterationCounter++;
      }
    }

    catch ( const FailedToConverge& e ) {
      throw ReturnMappingFailedException();
    }

    const MaterialState newMaterialState = extractMaterialState( X, trialState, activeSurfaces );

    return ReturnMapResult{ .materialState     = newMaterialState,
                            .dStrainPlastic    = computePlasticStrainIncrement( newMaterialState.stress,
                                                                             trialState.stress,
                                                                             elasticCompliance ),
                            .algorithmicModuli = computeAlgorithmicModuli( dR_dX_,
                                                                           elasticStiffness,
                                                                           elasticCompliance ) };
  }
  VectorXd UnistrandPlasticity::createResidualScaleVector( const VectorXd& LHS )
  {
    const auto nRows = LHS.rows();
    VectorXd   S( nRows );
    S.setOnes();
    return S;
  }

  UnistrandPlasticity::AlgorithmicModuli UnistrandPlasticity::computeAlgorithmicModuli(
    const MatrixXd& dF_dX,
    const Matrix6d& elasticStiffness,
    const Matrix6d& elasticCompliance )
  {

    const auto&       Cel = elasticStiffness;
    AlgorithmicModuli result;

    const auto      sizeEq = dF_dX.rows();
    Eigen::MatrixXd dYdE   = Eigen::MatrixXd::Identity( sizeEq, sizeEq );

    const Eigen::MatrixXd dXdE = dF_dX.colPivHouseholderQr().solve( dYdE );

    result.dStress_dStrain        = dXdE.block< nStress, nStress >( idxS, idxS ) * Cel;
    result.dStrainPlastic_dStrain = -elasticCompliance * result.dStress_dStrain + Matrix6d::Identity();

    return result;
  }
  Vector6d UnistrandPlasticity::computePlasticStrainIncrement( const Vector6d& stressNew,
                                                               const Vector6d& trialStress,
                                                               const Matrix6d& elasticCompliance )
  {
    return elasticCompliance * ( trialStress - stressNew );
  }

  UnistrandPlasticity::MaterialState UnistrandPlasticity::extractMaterialState( const Eigen::VectorXd&  X,
                                                                                const MaterialState&    trialState,
                                                                                const YieldSurfFlagArr& activeSurfaces )
  {
    // initialize new material state with trial state
    MaterialState newMaterialState = trialState;
    newMaterialState.stress        = X.segment< nStress >( idxS );

    int idxInternal_ = idxInternal;
    for ( int i = 0; i < nYieldSurfaces; i++ ) {
      // if ( activeSurfaces( i ) ) { 
      newMaterialState.alpha( i ) = X( idxInternal_ );
      idxInternal_ += 2;
     // }
    }

    return newMaterialState;
  };
  std::tuple< VectorXd, MatrixXd > UnistrandPlasticity::dR_dX( const VectorXd&                        X,
                                                               const MaterialState&                   trialState,
                                                               const Matrix6d&                        elasticStiffness,
                                                               UnistrandPlasticity::YieldSurfFlagArr& activeSurfaces )
  {

    using namespace Marmot::AutomaticDifferentiation;

    vector_to_vector_function_type_dual func = [&]( const autodiff::VectorXdual& X_ ) {
      return R( X_, trialState, elasticStiffness, activeSurfaces );
    };

    return AutomaticDifferentiation::dF_dX( func, X );
  }
} // namespace Marmot::Materials
