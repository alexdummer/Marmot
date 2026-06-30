#include "Marmot/GradientLinearElastic.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotPhaseFieldEnergyDegradation.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;

  GradientLinearElastic::GradientLinearElastic( const double* materialProperties,
                                                int           nMaterialProperties,
                                                int           materialNumber )
    : MarmotMaterialGradientPlasticityHypoElastic< 1 >( materialProperties, nMaterialProperties, materialNumber ),
      E( materialProperties[0] ),
      nu( materialProperties[1] ),
      lambda( E * nu / ( ( 1 + nu ) * ( 1 - 2 * nu ) ) ),
      mu( E / ( 2 * ( 1 + nu ) ) )
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
    using namespace Fastor;
    using namespace FastorStandardTensors;
    using namespace FastorIndices;
    const Tensor33d& _I = Spatial3D::I;
    res.stress += lambda * trace( inc.dStrain ) * _I + 2 * mu * inc.dStrain;
    res.f( 0 )          = 0;
    tan.dStressddStrain = lambda * einsum< ij, kl, to_ijkl >( _I, _I ) + 2 * mu * einsum< ik, jl, to_ijkl >( _I, _I );
  }

} // namespace Marmot::Materials
