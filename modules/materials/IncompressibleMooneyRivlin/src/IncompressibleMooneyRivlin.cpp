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
    using scalar  = autodiff::dual;
    const auto& F = deformation.F;

    using namespace ContinuumMechanics;

    // principal stretches lambda_i and eigenvectors Q of the right Cauchy-Green tensor C
    const auto [lambda, Q] = DeformationMeasures::principalStretches( F );
    const scalar J         = DeformationMeasures::volumeRatio( F );

    // the classical incompressible Mooney-Rivlin potential
    //   Psi_iso = C1*(I1bar-3) + C2*(I2bar-3)
    // is exactly the two-term Ogden model with (mu_1,alpha_1)=(2*C1,2), (mu_2,alpha_2)=(-2*C2,-2)
    // (using I2bar = sum_i lambda_bar_i^-2 for incompressible stretches) - reuse the validated
    // Ogden potential/stress rather than re-deriving the same math.
    const double C1 = materialProperties[0];
    const double C2 = materialProperties[1];

    const double mu[2]    = { 2.0 * C1, -2.0 * C2 };
    const double alpha[2] = { 2.0, -2.0 };

    const auto [psi_, tauPrincipal] = EnergyDensityFunctions::Incompressible::OgdenPrincipalKirchhoffStress( lambda,
                                                                                                             J,
                                                                                                             mu,
                                                                                                             alpha,
                                                                                                             2 );

    const auto tau = StressMeasures::KirchhoffStressFromPrincipalKirchhoffStress( tauPrincipal, lambda, Q, F );

    response.tau                  = tau;
    response.elasticEnergyDensity = Math::makeReal( psi_ );
    response.dissipation          = 0.0;
  }
} // namespace Marmot::Materials
