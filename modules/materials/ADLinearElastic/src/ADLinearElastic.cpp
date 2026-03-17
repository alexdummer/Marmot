#include "Marmot/ADLinearElastic.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotTypedefs.h"
#include "autodiff/forward/dual/eigen.hpp"

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace autodiff;
  using namespace Eigen;
  using namespace ContinuumMechanics::Elasticity;

  ADLinearElastic::ADLinearElastic( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialHypoElasticAD::MarmotMaterialHypoElasticAD( materialProperties,
                                                                nMaterialProperties,
                                                                materialNumber ),
      E( materialProperties[0] ),
      nu( materialProperties[1] )
  {
    assert( nMaterialProperties >= 2 );
  }

  double ADLinearElastic::getDensity()
  {
    if ( nMaterialProperties >= 3 )
      return materialProperties[2];
    throw std::runtime_error(
      std::string( MakeString() << __PRETTY_FUNCTION__ << ": Density not specified for ADLinearElastic." ) );
  }

  double ADLinearElastic::getDampingCoefficient()
  {
    if ( nMaterialProperties >= 4 )
      return materialProperties[3];
    throw std::runtime_error( std::string(
      MakeString() << __PRETTY_FUNCTION__ << ": Damping coefficient not specified for ADLinearElastic." ) );
  }
  void ADLinearElastic::computeStressAD( state3DAD&            state,
                                         const autodiff::dual* dStrain,
                                         const timeInfo&       timeInfo ) const
  {
    mVector6dual            s( state.stress );
    const mVector6dualConst dE( dStrain );

    const MatrixXdual C( ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu ) );

    s = s + C * dE;
  }
} // namespace Marmot::Materials
