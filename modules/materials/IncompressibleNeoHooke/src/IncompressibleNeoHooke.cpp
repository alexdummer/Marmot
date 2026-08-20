#include "Marmot/IncompressibleNeoHooke.h"
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

  IncompressibleNeoHooke::IncompressibleNeoHooke( const double* materialProperties,
                                                  int           nMaterialProperties,
                                                  int           materialLabel )
    : MarmotMaterialFiniteStrainAD( materialProperties, nMaterialProperties, materialLabel )
  {
    stateLayout.finalize();
  }

  void IncompressibleNeoHooke::computeStressAD( ConstitutiveResponseAD< 3 >& response,
                                                const DeformationAD< 3 >&    deformation,
                                                const TimeIncrement&         timeIncrement ) const
  {
    using scalar  = autodiff::dual;
    const auto& F = deformation.F;

    using namespace ContinuumMechanics;

    // principal stretches lambda_i and eigenvectors Q of the right Cauchy-Green tensor C
    const auto [lambda, Q] = DeformationMeasures::principalStretches( F );
    const scalar J         = DeformationMeasures::volumeRatio( F );

    // the classical incompressible Neo-Hookean potential is exactly the one-term Ogden model
    // Psi_iso = (mu/2) * (lambda_bar_1^2 + lambda_bar_2^2 + lambda_bar_3^2 - 3), i.e. mu=materialProperties[0],
    // alpha=2 - reuse the validated Ogden potential/stress rather than re-deriving the same math.
    const double mu[1]    = { materialProperties[0] };
    const double alpha[1] = { 2.0 };

    const auto [psi_, tauPrincipal] = EnergyDensityFunctions::Incompressible::OgdenPrincipalKirchhoffStress( lambda,
                                                                                                             J,
                                                                                                             mu,
                                                                                                             alpha,
                                                                                                             1 );

    const auto tau = StressMeasures::KirchhoffStressFromPrincipalKirchhoffStress( tauPrincipal, lambda, Q, F );

    response.tau                  = tau;
    response.elasticEnergyDensity = Math::makeReal( psi_ );
    response.dissipation          = 0.0;
  }
} // namespace Marmot::Materials
