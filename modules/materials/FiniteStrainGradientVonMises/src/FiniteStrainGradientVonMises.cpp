#include "Marmot/FiniteStrainGradientVonMises.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteStrainPlasticity.h"
#include "Marmot/MarmotStressMeasures.h"

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace FastorIndices;
  using namespace FastorStandardTensors;
  using namespace autodiff;
  using namespace ContinuumMechanics;

  FiniteStrainGradientVonMises::FiniteStrainGradientVonMises( const double* materialProperties,
                                                              int           nMaterialProperties,
                                                              int           materialNumber )
    : MarmotMaterialGradientPlasticityFiniteStrainAD< 1 >( materialProperties, nMaterialProperties, materialNumber ),
      K( materialProperties[0] ),
      G( materialProperties[1] ),
      fy0( materialProperties[2] ),
      H( materialProperties[3] ),
      g( materialProperties[4] )
  {
    stateLayout.add( "Fp", 9 );
    stateLayout.add( "kappa", 1 );
    stateLayout.add( "laplaceKappa", 1 );
    stateLayout.finalize();
  }

  double FiniteStrainGradientVonMises::getDensity( const double* stateVars ) const
  {
    if ( nMaterialProperties <= 5 )
      throw std::runtime_error( "Density not provided in material properties." );
    return materialProperties[5];
  }

  std::vector< double > FiniteStrainGradientVonMises::getNonlocalViscosity( const double* stateVars ) const
  {
    if ( nMaterialProperties <= 6 )
      throw std::runtime_error( "Viscosity not provided in material properties." );
    return { materialProperties[6] };
  }

  void FiniteStrainGradientVonMises::initializeYourself( double* stateVars, int nStateVars )
  {
    for ( int i = 0; i < nStateVars; ++i ) {
      stateVars[i] = 0.0;
    }

    TensorMap33d Fp = stateLayout.getAs< TensorMap33d >( stateVars, "Fp" );
    memcpy( Fp.data(), Spatial3D::I.data(), 9 * sizeof( double ) );
  }

  dual FiniteStrainGradientVonMises::fischerBurmeisterFunction( const dual a, const dual b, const double epsilon ) const
  {
    const dual sqrtTerm = sqrt( a * a + b * b + epsilon );
    const dual f        = sqrtTerm - ( a + b );
    return f;
  }

  dual FiniteStrainGradientVonMises::fy( const dual& kappa, const dual& laplaceKappa ) const
  {
    return fy0 + H * kappa - g * laplaceKappa;
  }

  void FiniteStrainGradientVonMises::computeStressAD( responseAD& res, const incrementAD& inc ) const
  {
    TensorMap33d      FpOld_map = stateLayout.getAs< TensorMap33d >( res.stateVars, "Fp" );
    Tensor33t< dual > FpOld     = Marmot::makeDual( Tensor33d( FpOld_map ) );

    const dual& dLambda = inc.dLambda( 0 );

    // Accumulate kappa and laplaceKappa from state variables
    double& kappaOld        = stateLayout.getAs< double& >( res.stateVars, "kappa" );
    double& laplaceKappaOld = stateLayout.getAs< double& >( res.stateVars, "laplaceKappa" );

    autodiff::dual kappa        = kappaOld + inc.dLambda( 0 );
    autodiff::dual laplaceKappa = laplaceKappaOld + inc.laplaceDLambda( 0 );

    // 1. Trial state: F^{e,trial} = F * (F^{p,old})^{-1}
    Tensor33t< dual > FeTrial = inc.deformation.F % Fastor::inverse( FpOld );
    Tensor33t< dual > CeTrial = transpose( FeTrial ) % FeTrial;

    auto [psiTrial, dPsi_dCeTrial] = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( CeTrial, K, G );
    Tensor33t< dual > STrial       = dPsi_dCeTrial + dPsi_dCeTrial; //  2 * dPsi_dCeTrial

    dual              J        = determinant( inc.deformation.F );
    dual              Jinv     = dual( 1.0 ) / J;
    Tensor33t< dual > tauTrial = FeTrial % STrial % transpose( FeTrial );
    Tensor33t< dual > tTrial   = multiplyFastorTensorWithScalar( tauTrial, Jinv );

    // Deviatoric trial Cauchy stress and J2
    Tensor33t< dual > tTrial_dev = Marmot::deviatoric( tTrial );
    dual              J2_tTrial  = dual( 0.5 ) * einsum_ij_ij_hardcoded( tTrial_dev, tTrial_dev );
    dual              f_tr       = sqrt( dual( 3.0 ) * J2_tTrial ) - fy( kappa, laplaceKappa );

    // Flow direction: df_p/dt = (3/2) * t_dev / sqrt(3 J2)
    // df_p/dM_{KL} = (df_p/dt_{ij}) * (1/J) * F^{e,-1}_{Ki} * F^{e}_{jL}
    Tensor33t< dual > dfp_dt;
    if ( Math::makeReal( J2_tTrial ) > 1e-12 ) {
      dual factor = dual( 1.5 ) / sqrt( dual( 3.0 ) * J2_tTrial );
      dfp_dt      = multiplyFastorTensorWithScalar( tTrial_dev, factor );
    }
    else {
      dfp_dt = makeDual( Tensor33d( 0.0 ) );
    }

    Tensor33t< dual > FeTrial_inv     = Fastor::inverse( FeTrial );
    Tensor33t< dual > dfp_dM_unscaled = Tensor33t< dual >( transpose( FeTrial_inv ) % dfp_dt % transpose( FeTrial ) );
    Tensor33t< dual > dfp_dM          = multiplyFastorTensorWithScalar( dfp_dM_unscaled, Jinv );

    // 2. Exponential map: delta Fp = exp( dLambda * dfp_dM )
    Tensor33t< dual > dGp     = multiplyFastorTensorWithScalar( dfp_dM, dLambda );
    Tensor33t< dual > deltaFp = ContinuumMechanics::FiniteStrain::Plasticity::FlowIntegration::exponentialMap( dGp );

    Tensor33t< dual > FpNew = deltaFp % FpOld;

    // 3. Updated elastic state: F^e = F * (F^{p,new})^{-1}
    Tensor33t< dual > Fe  = inc.deformation.F % Fastor::inverse( FpNew );
    Tensor33t< dual > Ce  = transpose( Fe ) % Fe;
    auto [psi, dPsi_dCe]  = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( Ce, K, G );
    Tensor33t< dual > S   = multiplyFastorTensorWithScalar( dPsi_dCe, dual( 2.0 ) );
    Tensor33t< dual > tau = Fe % S % transpose( Fe );

    // 4. Fischer-Burmeister complementarity function
    double scale = 1e6;
    dual   fFB   = fischerBurmeisterFunction( -f_tr, dLambda * scale, 1e-12 );

    // 5. Update state variables (primal values only)
    memcpy( stateLayout.getPtr( res.stateVars, "Fp" ), makeReal( FpNew ).data(), 9 * sizeof( double ) );
    stateLayout.getAs< double& >( res.stateVars, "kappa" )        = Math::makeReal( kappa );
    stateLayout.getAs< double& >( res.stateVars, "laplaceKappa" ) = Math::makeReal( laplaceKappa );

    // 6. Write response
    res.tau                  = tau;
    res.f( 0 )               = fFB;
    res.elasticEnergyDensity = psi;
    res.dissipation          = dual( 0.0 );
  }
} // namespace Marmot::Materials
