#include "Marmot/GradientVonMises.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotTypedefs.h"
#include <algorithm>
#include <limits>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;

  GradientVonMises::GradientVonMises( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialGradientPlasticityHypoElastic( materialProperties, nMaterialProperties, materialNumber )
  {
    stateLayout.add( "kappaLocal", 1 );
    stateLayout.finalize();
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
      throw std::runtime_error( "Nonlocal viscosity not provided in material properties for GradientVonMises." );
    return { materialProperties[6] };
  }

  void GradientVonMises::computeStress( response& res, tangents& tan, const increment& inc ) const
  {
    const double& E         = materialProperties[0];
    const double& nu        = materialProperties[1];
    const double& sigmaY0   = materialProperties[2];
    const double& HNonlocal = materialProperties[3];
    const double& l         = materialProperties[4];

    const Matrix6d Cel = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
    const double   G   = E / ( 2. * ( 1. + nu ) );

    double& kappaLocal = stateLayout.getAs< double& >( res.stateVars, "kappaLocal" );

    const Vector6d trialStress = res.stress + Cel * inc.dStrain;
    using namespace ContinuumMechanics::VoigtNotation;
    const double rhoTrial = std::sqrt( 2. * Invariants::J2( trialStress ) );

    const double KNonlocal = inc.K( 0 );
    const double fTrial    = rhoTrial - Constants::sqrt2_3 * ( sigmaY0 + HNonlocal * KNonlocal );

    res.c( 0 )      = l * l;
    res.KLocal( 0 ) = kappaLocal;

    tan.dStressddStrain = Cel;
    tan.dStressddK.setZero();
    tan.dKLocalddStrain.setZero();
    tan.dKLocalddK.setZero();
    tan.dcddK.setZero();
    tan.d2cddK2.setZero();

    if ( fTrial <= 0. || rhoTrial < 1e-14 ) {
      res.stress = trialStress;
      return;
    }

    const Vector6d n       = IDev * trialStress / rhoTrial;
    const double   dKappa  = std::max( 0.0,
                                    fTrial /
                                      std::max( Constants::sqrt6 * G, std::numeric_limits< double >::epsilon() ) );
    const double   dLambda = Constants::sqrt3_2 * dKappa;

    res.stress = trialStress - 2. * G * dLambda * n;

    kappaLocal += dKappa;
    res.KLocal( 0 ) = kappaLocal;

    Matrix6d IDevHalfShear = IDev;
    IDevHalfShear.block< 6, 3 >( 0, 3 ) *= 0.5;

    tan.dStressddStrain     = Cel - 2. * G * ( n * n.transpose() ) - 4. * G * G * dLambda / rhoTrial * IDevHalfShear;
    tan.dStressddK.col( 0 ) = Constants::sqrt2_3 * HNonlocal * n;

    tan.dKLocalddStrain.row( 0 ) = Constants::sqrt2_3 * n.transpose();
    tan.dKLocalddK( 0, 0 )       = -HNonlocal / ( 3. * G );
  }

} // namespace Marmot::Materials
