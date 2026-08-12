#include "Marmot/GradientVonMises.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotPhaseFieldEnergyDegradation.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;

  GradientVonMises::GradientVonMises( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialGradientPlasticityHypoElastic< 1 >( materialProperties, nMaterialProperties, materialNumber ),
      E( materialProperties[0] ),
      C( ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( materialProperties[0], materialProperties[1] ) ),
      fy0( materialProperties[2] ),
      H( materialProperties[3] ),
      g( materialProperties[4] ),
      implementation( materialProperties[5] == 0 ? Implementation::standard : Implementation::fischer_burmeister )
  {
    // optional, so that existing property vectors of six to eight entries keep working
    if ( nMaterialProperties > 8 && materialProperties[8] > 0.0 )
      fbSmoothingRelative = materialProperties[8];

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
    // map response and increment variables for easier access
    mVector6d       stress( res.stress.data() );
    const Vector6d& dStrain = inc.dStrain;
    double&         f       = res.f( 0 ); // yield function value

    // map to tangents
    mMatrix6d                     dStressddStrain( tan.dStressddStrain.data() );
    Map< Matrix< double, 6, 1 > > dStressddLambda( tan.dStressddLambda.data() );
    Map< Matrix< double, 1, 6 > > dF_ddStrain( tan.dFddStrain.data() );
    double&                       dF_dKappa        = tan.dFddLambda( 0, 0 );
    double&                       dF_dLaplaceKappa = tan.dFddLaplacian( 0, 0 );

    // get state variables
    double& kappa        = stateLayout.getAs< double& >( res.stateVars, "kappa" );
    double& laplaceKappa = stateLayout.getAs< double& >( res.stateVars, "laplaceKappa" );

    // update kappa and laplaceKappa
    kappa += inc.dLambda( 0 );
    laplaceKappa += inc.laplaceDLambda( 0 );

    // compute trial stress
    Vector6d trialStress = stress + C * dStrain;
    // handle zero increment
    if ( inc.dStrain.isZero( 1e-14 ) && inc.laplaceDLambda.isZero( 1e-14 ) && inc.dLambda.isZero( 1e-14 ) ) {
      stress           = trialStress;
      f                = 0; // ensure yield function is exactly zero for zero increment
      dF_dKappa        = E; // set hardening derivative to Young's modulus for zero increment
      dF_dLaplaceKappa = g; // set gradient hardening derivative to zero for zero increment
      dF_ddStrain.setZero();
      dStressddStrain = C;
      // std::cout << "Zero increment detected, skipping return mapping and setting stress to trial stress." <<
      // std::endl; std::exit( 0 ); //
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
      dF_ddStrain      = dF_dStress_tr.transpose() * C;
      dStressddStrain  = C;
      return;
    }
    // update stress with trial return mapping direction
    stress          = trialStress - C * ( dLambda * dF_dStress_tr );
    dStressddStrain = C - dLambda * C * d2F_dStress2_tr * C;
    dStressddLambda = -C * dF_dStress_tr;
    Vector6d dF_dStress;
    Matrix6d d2F_dStress2;
    std::tie( f, dF_dStress, d2F_dStress2, dF_dKappa, dF_dLaplaceKappa ) = yieldFunction( stress, kappa, laplaceKappa );
    dF_ddStrain                                                          = dF_dStress.transpose() * dStressddStrain;
    dF_dKappa += dF_dStress.dot( dStressddLambda );
  }

  void GradientVonMises::computeStressFischerBurmeister( response& res, tangents& tan, const increment& inc ) const
  {

    // map response and increment variables for easier access
    mVector6d       stress( res.stress.data() );
    const Vector6d& dStrain = inc.dStrain;
    double&         f       = res.f( 0 ); // yield function value

    // map to tangents
    mMatrix6d                     dStressddStrain( tan.dStressddStrain.data() );
    Map< Matrix< double, 6, 1 > > dStressddLambda( tan.dStressddLambda.data() );
    Map< Matrix< double, 1, 6 > > dF_ddStrain( tan.dFddStrain.data() );
    double&                       dF_dKappa        = tan.dFddLambda( 0, 0 );
    double&                       dF_dLaplaceKappa = tan.dFddLaplacian( 0, 0 );

    // get state variables
    double& kappa        = stateLayout.getAs< double& >( res.stateVars, "kappa" );
    double& laplaceKappa = stateLayout.getAs< double& >( res.stateVars, "laplaceKappa" );

    // update kappa and laplaceKappa
    kappa += inc.dLambda( 0 );
    laplaceKappa += inc.laplaceDLambda( 0 );

    // compute trial stress
    Vector6d trialStress = stress + C * dStrain;
    // handle zero increment
    // if ( inc.dStrain.norm() == 0 && inc.laplaceDLambda( 0 ) == 0 && inc.dLambda( 0 ) == 0 ) {
    //   stress = trialStress;
    //   Vector6d dF_dStress;
    //   Matrix6d d2F_dStress2;
    //   std::tie( f, dF_dStress, d2F_dStress2, dF_dKappa, dF_dLaplaceKappa ) = yieldFunction( stress,
    //                                                                                         kappa,
    //                                                                                         laplaceKappa );
    //   dF_ddStrain                                                          = dF_dStress.transpose() * C;
    //   dStressddStrain                                                      = C;
    //   return;
    // }
    const double& dLambda = inc.dLambda( 0 );

    auto [f_tr, dF_dStress_tr, d2F_dStress2_tr, dF_dKappa_tr, dF_dLaplaceKappa_tr] = yieldFunction( trialStress,
                                                                                                    kappa,
                                                                                                    laplaceKappa );
    // update stress with trial return mapping direction
    stress = trialStress - C * ( dLambda * dF_dStress_tr );

    dStressddStrain = C - dLambda * C * d2F_dStress2_tr * C;
    dStressddLambda = -C * dF_dStress_tr;

    Vector6d dF_dStress;
    Matrix6d d2F_dStress2;
    std::tie( f_tr, dF_dStress_tr, d2F_dStress2_tr, dF_dKappa_tr, dF_dLaplaceKappa_tr ) = yieldFunction( stress,
                                                                                                         kappa,
                                                                                                         laplaceKappa );

    double df_da, df_db;
    double scale = 1e4; // scaling factor to improve conditioning of the Fischer-Burmeister function derivatives

    // The corner of the complementarity function has to be blunt enough that a residual sized
    // perturbation of the yield function does not swing its derivative, see fbSmoothingRelative.
    // Scaled by the yield strength because the first argument is a stress.
    const double epsilon = std::pow( fbSmoothingRelative * fy0, 2 );

    std::tie( f,
              df_da,
              df_db ) = fischerBurmeisterFunction( -f_tr,
                                                   dLambda * scale,
                                                   epsilon ); // using Fischer-Burmeister to enforce yield condition

                                                              // Compute derivatives of the Fischer-Burmeister function
    dF_dStress = -df_da * dF_dStress_tr;
    dF_dKappa  = -df_da * dF_dKappa_tr + df_db * scale;

    dF_ddStrain = dF_dStress.transpose() * dStressddStrain;
    dF_dKappa += dF_dStress.dot( dStressddLambda );
    dF_dLaplaceKappa = -df_da * dF_dLaplaceKappa_tr + df_db * 0;
  }

  std::tuple< double, double, double > GradientVonMises::fy( double kappa, double laplaceKappa ) const
  {
    const double sigmaY = fy0 + H * kappa - g * laplaceKappa; // yield stress as a function of kappa
    return { sigmaY, H, -g };
  }

  // compute the von Mises yield function value and its derivatives with respect to stress, kappa, and laplaceKappa
  std::tuple< double, Vector6d, Matrix6d, double, double > GradientVonMises::yieldFunction(
    const Vector6d& stress,
    const double&   kappa,
    const double&   laplaceKappa ) const
  {
    using namespace Marmot::ContinuumMechanics::VoigtNotation;
    const auto [sigmaY, dSigmaY_dKappa, dSigmaY_dLaplaceKappa] = fy( kappa, laplaceKappa );
    const double J2                                            = Invariants::J2( stress );
    const double f                                             = std::sqrt( 3.0 * J2 ) - sigmaY; // yield

    if ( J2 < 1e-12 ) {
      return { f, Vector6d::Zero(), Matrix6d::Zero(), -dSigmaY_dKappa, -dSigmaY_dLaplaceKappa };
    }

    // 3. Compute base derivatives if J2 is safely non-zero
    const Vector6d dJ2_dStress   = Derivatives::dJ2_dStress( stress );
    const Matrix6d d2J2_dStress2 = Derivatives::d2J2_dStress2( stress );

    // 4. Compute correct yield function derivatives
    const double   root3J2    = std::sqrt( 3.0 * J2 );
    const Vector6d dF_dStress = dJ2_dStress * ( 3.0 / ( 2.0 * root3J2 ) );

    // Corrected the math scalar multiplier here from 9/8 to sqrt(3)/4
    const double   scalar2Matrix = std::sqrt( 3.0 ) / ( 4.0 * std::pow( J2, 1.5 ) );
    const Matrix6d d2F_dStress2  = d2J2_dStress2 * ( 3.0 / ( 2.0 * root3J2 ) ) -
                                  ( dJ2_dStress * dJ2_dStress.transpose() ) * scalar2Matrix;

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
