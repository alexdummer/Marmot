#include "Marmot/GradientLinearElastic.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotPhaseFieldEnergyDegradation.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;

  GradientLinearElastic::GradientLinearElastic( const double* materialProperties,
                                                int           nMaterialProperties,
                                                int           materialNumber )
    : MarmotMaterialGradientPlasticityHypoElastic< 1 >( materialProperties, nMaterialProperties, materialNumber ),
      C( ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( materialProperties[0], materialProperties[1] ) )
  {
    initializeStateLayout();
  }

  double GradientLinearElastic::getDensity( const double* stateVars ) const
  {
    if ( nMaterialProperties <= 2 )
      throw std::runtime_error( "Density not provided in material properties for GradientLinearElastic." );
    return materialProperties[2];
  }

  std::vector< double > GradientLinearElastic::getNonlocalViscosity( const double* stateVars ) const
  {
    if ( nMaterialProperties <= 3 )
      throw std::runtime_error( "Viscosity not provided in material properties for GradientLinearElastic." );
    return { materialProperties[3] };
  }
  void GradientLinearElastic::computeStress( response& res, tangents& tan, const increment& inc ) const
  {
    res.stress += C * inc.dStrain;
    res.f.setZero();
    tan.dStressddStrain = C;
    tan.dStressddLambda.setZero();
    tan.dStressddLaplacian.setZero();
  }

} // namespace Marmot::Materials
