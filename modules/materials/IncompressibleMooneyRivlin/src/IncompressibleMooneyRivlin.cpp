#include "Marmot/IncompressibleMooneyRivlin.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStressMeasures.h"
#include <Fastor/expressions/linalg_ops/unary_trans_op.h>
#include <Fastor/tensor/Tensor.h>
#include <Fastor/tensor_algebra/einsum_explicit.h>
#include <Fastor/tensor_algebra/indicial.h>
#include <autodiff/forward/dual/dual.hpp>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace FastorIndices;
  using namespace FastorStandardTensors;

  IncompressibleMooneyRivlin::IncompressibleMooneyRivlin( const double* materialProperties,
                                                          int           nMaterialProperties,
                                                          int           materialLabel )
    : MarmotMaterialFiniteStrainAD( materialProperties, nMaterialProperties, materialLabel )
  {
    stateLayout.finalize();
  }

  void IncompressibleMooneyRivlin::computeStressAD( ConstitutiveResponseAD< 3 >& response,
                                                    const DeformationAD< 3 >&    deformation,
                                                    const TimeIncrement&         timeIncrement ) const
  {
    const double& C1 = materialProperties[0];
    const double& C2 = materialProperties[1];
    const auto&   F_ = deformation.F;

    using namespace ContinuumMechanics;

    // Unlike Ogden's general (non-integer) exponents, the classical Mooney-Rivlin potential is
    // quadratic in the isochoric stretches and can therefore be expressed directly in terms of the
    // invariants of C - no spectral decomposition of C is required.
    const auto C = DeformationMeasures::rightCauchyGreen( F_ );

    const auto [psi_, dPsi_dC] = EnergyDensityFunctions::Incompressible::FirstOrderDerived::MooneyRivlinPotential( C,
                                                                                                                   C1,
                                                                                                                   C2 );

    Tensor33t< autodiff::dual > PK2 = multiplyFastorTensorWithScalar( dPsi_dC, autodiff::dual( 2.0 ) );

    const auto tau                = StressMeasures::KirchhoffStressFromPK2( PK2, F_ );
    response.tau                  = tau;
    response.elasticEnergyDensity = Math::makeReal( psi_ );
    response.dissipation          = 0.0;
  }
} // namespace Marmot::Materials
