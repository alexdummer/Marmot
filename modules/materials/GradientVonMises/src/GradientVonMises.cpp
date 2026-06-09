#include "Marmot/GradientVonMises.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotPhaseFieldEnergyDegradation.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;

  GradientVonMises::GradientVonMises( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialGradientPlasticityHypoElastic< 1 >( materialProperties, nMaterialProperties, materialNumber ),
      C( ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( materialProperties[0], materialProperties[1] ) ),
      fy0( materialProperties[2] ),
      H( materialProperties[3] ),
      g( materialProperties[4] )
  {
    initializeStateLayout();
  }

  double GradientVonMises::getDensity( const double* stateVars ) const
  {
    if ( nMaterialProperties <= 5 )
      throw std::runtime_error( "Density not provided in material properties for GradientVonMises." );
    return materialProperties[5];
  }

  std::vector< double > GradientVonMises::getNonlocalViscosity( const double* stateVars ) const
  {
    if ( nMaterialProperties <= 6 )
      throw std::runtime_error( "Viscosity not provided in material properties for GradientVonMises." );
    return { materialProperties[6] };
  }
  void GradientVonMises::computeStress( response& res, tangents& tan, const increment& inc ) const
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
    if ( inc.dStrain.norm() == 0 && inc.laplaceDLambda( 0 ) == 0 && inc.dLambda( 0 ) == 0 ) {
      // elastic step
      stress = trialStress;
      Vector6d dF_dStress;
      Matrix6d d2F_dStress2;
      std::tie( f, dF_dStress, d2F_dStress2, dF_dKappa, dF_dLaplaceKappa ) = yieldFunction( stress,
                                                                                            kappa,
                                                                                            laplaceKappa );
      dF_ddStrain                                                          = dF_dStress.transpose() * C;
      dStressddStrain                                                      = C;
      return;
    }
    const double& dLambda = inc.dLambda( 0 );

    const auto [f_tr,
                dF_dStress_tr,
                d2F_dStress2_tr,
                dF_dKappa_tr,
                dF_dLaplaceKappa_tr] = yieldFunction( trialStress, kappa, laplaceKappa );
    // update stress with trial return mapping direction
    stress = trialStress - C * ( dLambda * dF_dStress_tr );

    dStressddStrain = C - dLambda * C * d2F_dStress2_tr * C;
    dStressddLambda = -C * dF_dStress_tr;

    Vector6d dF_dStress;
    Matrix6d d2F_dStress2;
    std::tie( f, dF_dStress, d2F_dStress2, dF_dKappa, dF_dLaplaceKappa ) = yieldFunction( stress, kappa, laplaceKappa );
    dF_ddStrain                                                          = dF_dStress.transpose() * dStressddStrain;
    dF_dKappa += dF_dStress.transpose() * dStressddLambda;
  }

} // namespace Marmot::Materials
