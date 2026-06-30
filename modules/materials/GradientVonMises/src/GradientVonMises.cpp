#include "Marmot/GradientVonMises.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotPhaseFieldEnergyDegradation.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Fastor;
  using namespace Marmot;
  using namespace Marmot::FastorStandardTensors;
  using namespace Marmot::FastorIndices;

  GradientVonMises::GradientVonMises( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialGradientPlasticityHypoElastic< 1 >( materialProperties, nMaterialProperties, materialNumber ),
      E( materialProperties[0] ),
      nu( materialProperties[1] ),
      lambda( E * nu / ( ( 1. - 2. * nu ) * ( 1. + nu ) ) ),
      mu( E / ( 2. * ( 1. + nu ) ) ),
      fy0( materialProperties[2] ),
      H( materialProperties[3] ),
      g( materialProperties[4] ),
      implementation( materialProperties[5] == 0 ? Implementation::standard : Implementation::fischer_burmeister )
  {
    initializeStateLayout();
  }

  double GradientVonMises::getDensity( const double* stateVars ) const
  {
    if ( nMaterialProperties <= 6 )
      throw std::runtime_error( "Density not provided in material properties for GradientVonMises." );
    return materialProperties[6];
  }

  std::vector< double > GradientVonMises::getNonlocalViscosity( const double* stateVars ) const
  {
    if ( nMaterialProperties <= 7 )
      throw std::runtime_error( "Viscosity not provided in material properties for GradientVonMises." );
    return { materialProperties[7] };
  }
  void GradientVonMises::computeStress( response& res, tangents& tan, const increment& inc ) const
  {
    switch ( implementation ) {
    case Implementation::standard: computeStressStandard( res, tan, inc ); break;
    case Implementation::fischer_burmeister: computeStressFischerBurmeister( res, tan, inc ); break;
    default: throw std::runtime_error( "Invalid implementation choice for GradientVonMises." );
    }
  }

  void GradientVonMises::computeStressStandard( response& res, tangents& tan, const increment& inc ) const
  {
    using namespace Fastor;
    using namespace Marmot::FastorStandardTensors;
    // map response and increment variables for easier access
    TensorMap33d stress( res.stress.data() );
    double&      f = res.f( 0 ); // yield function value

    // map to tangents
    TensorMap3333d dStressddStrain( tan.dStressddStrain.data() );
    TensorMap33d   dStressddLambda( tan.dStressddLambda.data() );
    TensorMap33d   dF_ddStrain( tan.dFddStrain.data() );
    double&        dF_dKappa        = tan.dFddLambda( 0, 0 );
    double&        dF_dLaplaceKappa = tan.dFddLaplacian( 0, 0 );

    // get state variables
    double& kappa        = stateLayout.getAs< double& >( res.stateVars, "kappa" );
    double& laplaceKappa = stateLayout.getAs< double& >( res.stateVars, "laplaceKappa" );

    // update kappa and laplaceKappa
    kappa += inc.dLambda( 0 );
    laplaceKappa += inc.laplaceDLambda( 0 );

    // compute trial stress
    Tensor33d trialStress;
    trialStress     = stress + lambda * trace( inc.dStrain ) * Spatial3D::I + 2 * mu * inc.dStrain;
    Tensor3333d Cel = lambda * outer( Spatial3D::I, Spatial3D::I ) + 2 * mu * Spatial3D::I4;
    // handle zero increment
    if ( norm( inc.dStrain ) < 1e-14 && norm( inc.laplaceDLambda ) < 1e-14 && norm( inc.dLambda ) < 1e-14 ) {
      stress           = trialStress;
      f                = 0; // ensure yield function is exactly zero for zero increment
      dF_dKappa        = E; // set hardening derivative to Young's modulus for zero increment
      dF_dLaplaceKappa = g; // set gradient hardening derivative to zero for zero increment
      dF_ddStrain      = Tensor33d( 0. );
      dStressddStrain  = Cel;
      return;
    }
    const double& dLambda = inc.dLambda( 0 );

    const auto [f_tr,
                dF_dStress_tr,
                d2F_dStress2_tr,
                dF_dKappa_tr,
                dF_dLaplaceKappa_tr] = yieldFunction( trialStress, kappa, laplaceKappa );
    // elastic step
    if ( f_tr <= 0 ) {
      stress           = trialStress;
      res.f( 0 )       = 0; // ensure yield function is exactly zero for elastic step
      dF_dKappa        = E; // set hardening derivative to Young's modulus for elastic step
      dF_dLaplaceKappa = dF_dLaplaceKappa_tr;
      dF_ddStrain      = einsum< ij, ijkl >( dF_dStress_tr, Cel );
      dStressddStrain  = Cel;
      return;
    }
    // update stress with trial return mapping direction
    stress          = trialStress - Cel * ( dLambda * dF_dStress_tr );
    dStressddStrain = Cel - dLambda * Cel * d2F_dStress2_tr * Cel;
    dStressddLambda = -Cel * dF_dStress_tr;
    Tensor33d   dF_dStress;
    Tensor3333d d2F_dStress2;
    std::tie( f, dF_dStress, d2F_dStress2, dF_dKappa, dF_dLaplaceKappa ) = yieldFunction( stress, kappa, laplaceKappa );
    dF_ddStrain = einsum< ij, ijkl >( dF_dStress, dStressddStrain );
    dF_dKappa += einsum< ij, ij >( dF_dStress, dStressddLambda ).toscalar();
  }

  void GradientVonMises::computeStressFischerBurmeister( response& res, tangents& tan, const increment& inc ) const
  {

    // map response and increment variables for easier access
    TensorMap33d stress( res.stress.data() );
    double&      f = res.f( 0 ); // yield function value

    // map to tangents
    TensorMap3333d dStressddStrain( tan.dStressddStrain.data() );
    TensorMap33d   dStressddLambda( tan.dStressddLambda.data() );
    TensorMap33d   dF_ddStrain( tan.dFddStrain.data() );
    double&        dF_dKappa        = tan.dFddLambda( 0, 0 );
    double&        dF_dLaplaceKappa = tan.dFddLaplacian( 0, 0 );

    // get state variables
    double& kappa        = stateLayout.getAs< double& >( res.stateVars, "kappa" );
    double& laplaceKappa = stateLayout.getAs< double& >( res.stateVars, "laplaceKappa" );

    // update kappa and laplaceKappa
    kappa += inc.dLambda( 0 );
    laplaceKappa += inc.laplaceDLambda( 0 );

    // compute trial stress
    Tensor33d     trialStress = stress + lambda * trace( inc.dStrain ) * Spatial3D::I + 2 * mu * inc.dStrain;
    Tensor3333d   Cel         = lambda * outer( Spatial3D::I, Spatial3D::I ) + 2 * mu * Spatial3D::I4;
    const double& dLambda     = inc.dLambda( 0 );

    auto [f_tr, dF_dStress_tr, d2F_dStress2_tr, dF_dKappa_tr, dF_dLaplaceKappa_tr] = yieldFunction( trialStress,
                                                                                                    kappa,
                                                                                                    laplaceKappa );
    // update stress with trial return mapping direction
    stress = trialStress - Cel * ( dLambda * dF_dStress_tr );

    dStressddStrain = Cel - dLambda * Cel * d2F_dStress2_tr * Cel;
    dStressddLambda = -Cel * dF_dStress_tr;

    Tensor33d dF_dStress;
    // Tensor3333d d2F_dStress2;
    std::tie( f_tr, dF_dStress_tr, d2F_dStress2_tr, dF_dKappa_tr, dF_dLaplaceKappa_tr ) = yieldFunction( stress,
                                                                                                         kappa,
                                                                                                         laplaceKappa );

    double df_da, df_db;
    double scale      = 1e4; // scaling factor to improve conditioning of the Fischer-Burmeister function derivatives
    std::tie( f,
              df_da,
              df_db ) = fischerBurmeisterFunction( -f_tr,
                                                   dLambda * scale,
                                                   1e-16 ); // using Fischer-Burmeister to enforce yield condition

    // Compute derivatives of the Fischer-Burmeister function
    dF_dStress = -df_da * dF_dStress_tr;
    dF_dKappa  = -df_da * dF_dKappa_tr + df_db * scale;

    dF_ddStrain = einsum< ij, ijkl >( dF_dStress, dStressddStrain );
    dF_dKappa += einsum< ij, ij >( dF_dStress, dStressddLambda ).toscalar();
    dF_dLaplaceKappa = -df_da * dF_dLaplaceKappa_tr + df_db * 0;
  }

  std::tuple< double, double, double > GradientVonMises::fy( double kappa, double laplaceKappa ) const
  {
    const double sigmaY = fy0 + H * kappa - g * laplaceKappa; // yield stress as a function of kappa
    return { sigmaY, H, -g };
  }

  // compute the von Mises yield function value and its derivatives with respect to stress, kappa, and laplaceKappa
  std::tuple< double, Tensor33d, Tensor3333d, double, double > GradientVonMises::yieldFunction(
    const Tensor33d& stress,
    const double&    kappa,
    const double&    laplaceKappa ) const
  {

    const Tensor33d& _I = Spatial3D::I;

    const auto [sigmaY, dSigmaY_dKappa, dSigmaY_dLaplaceKappa] = fy( kappa, laplaceKappa );
    const Tensor33d devStress                                  = deviatoric( stress );
    const double    sigmaV = Marmot::Constants::sqrt3_2 * std::sqrt( einsum_ij_ij_hardcoded( devStress, devStress ) );
    double          f      = sigmaV - sigmaY; // yield

    if ( sigmaV < 1e-12 ) {
      return { f, 0.0, 0.0, -dSigmaY_dKappa, -dSigmaY_dLaplaceKappa };
    }

    // 3. Compute base derivatives if J2 is safely non-zero
    const Tensor33d   dF_dStress   = Marmot::Constants::sqrt3_2 * devStress / sigmaV;
    const Tensor3333d d2F_dStress2 = Marmot::Constants::sqrt3_2 *
                                     ( -1.0 / std::pow( sigmaV, 3 ) * outer( devStress, devStress ) +
                                       1. / sigmaV *
                                         ( einsum< ik, jl >( _I, _I ) - 1. / 3 * einsum< ij, IK, IL >( _I, _I, _I ) ) );

    return { f,
             dF_dStress,
             d2F_dStress2,
             -dSigmaY_dKappa,
             -dSigmaY_dLaplaceKappa }; // return yield function value and its derivatives
  }

  std::tuple< double, double, double > GradientVonMises::fischerBurmeisterFunction( const double a,
                                                                                    const double b,
                                                                                    const double epsilon ) const
  {
    const double sqrtTerm = std::sqrt( a * a + b * b + epsilon );
    const double f        = sqrtTerm - ( a + b );
    const double df_da    = a / sqrtTerm - 1.0;
    const double df_db    = b / sqrtTerm - 1.0;
    return { f, df_da, df_db };
  }

} // namespace Marmot::Materials
