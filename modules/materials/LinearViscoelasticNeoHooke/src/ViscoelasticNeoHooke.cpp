#include "Marmot/ViscoelasticNeoHooke.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotStressMeasures.h"
#include <Fastor/expressions/linalg_ops/unary_trans_op.h>
#include <Fastor/tensor/Tensor.h>
#include <Fastor/tensor_algebra/einsum_explicit.h>
#include <Fastor/tensor_algebra/indicial.h>
#include <autodiff/forward/dual/dual.hpp>
#include <map>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace FastorIndices;
  using namespace FastorStandardTensors;

  ViscoelasticNeoHooke::ViscoelasticNeoHooke( const double* materialProperties,
                                              int           nMaterialProperties,
                                              int           materialLabel )
    : MarmotMaterialFiniteStrain( materialProperties, nMaterialProperties, materialLabel ),
      K( materialProperties[0] ),
      G( materialProperties[1] ),
      nMaxwellElements( ( nMaterialProperties - 2 ) / 2 )
  {
    if ( nMaterialProperties < 4 )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                                << ": Not sufficient "
                                                   "materialProperties!" );

    // check if nMaterialProperties is divisible by 2
    if ( ( nMaterialProperties - 2 ) % 2 != 0 ) {
      throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": Wrong number of material properties!" );
    }

    // read in gamma and tau pairs
    tau   = Eigen::VectorXd( nMaxwellElements );
    gamma = Eigen::VectorXd( nMaxwellElements );
    for ( size_t i = 0; i < nMaxwellElements; ++i ) {
      gamma[i] = materialProperties[2 + i * 2];
      tau[i]   = materialProperties[2 + i * 2 + 1];
    }
  }

  void ViscoelasticNeoHooke::computeStress( ConstitutiveResponse< 3 >& response,
                                            AlgorithmicModuli< 3 >&    tangents,
                                            const Deformation< 3 >&    deformation,
                                            const TimeIncrement&       timeIncrement )
  {
    // deformation gradient
    const auto& F_ = deformation.F;

    using namespace ContinuumMechanics;
    // compute Cauchy-Green deformation
    const auto [C, dC_dF] = DeformationMeasures::FirstOrderDerived::CauchyGreen( F_ );

    // compute energy density, first and second partial derivatives wrt Cauchy Green deformation
    const auto [psi_, dPsi_dC, d2Psi_dCdC] = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( C, K, G );

    // compute PK2 of the purely elastic part
    Tensor33d   PK2     = 2. * dPsi_dC;
    Tensor3333d dPk2_dC = 2. * d2Psi_dCdC;

    // extract previous PK2 stress in the elastic part
    auto& PK2_old = stateVars->S0;

    // increment of PK2 stress in the elastic part
    Tensor33d dPK2 = PK2 - PK2_old;
    PK2_old        = PK2;

    const double& dt = timeIncrement.dT;
    // add contributions from Maxwell elements
    for ( size_t i = 0; i < nMaxwellElements; ++i ) {

      // extract state variables of the i-th Maxwell element
      Fastor::TensorMap< double, 3, 3 > h( stateVars->maxwellStateVars + i * 9 );

      // exponential decay factor
      const double expFactor = std::exp( -dt / tau[i] );

      // update state variable
      h = expFactor * h + ( 1.0 - expFactor ) * gamma[i] * tau[i] * dPK2 / dt;

      // add contribution to PK2 stress
      PK2 += h;

      // add contribution to tangent
      dPk2_dC += ( 1.0 - expFactor ) * gamma[i] * dPk2_dC;
    }

    // compute Kirchhoff stress
    const auto [kirchhoff, dTau_dPK2, dTau_dF] = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, F_ );
    response.tau                               = kirchhoff;
    response.rho                               = 1.0;
    response.elasticEnergyDensity              = psi_;

    // compute tangent operator
    tangents.dTau_dF = einsum< ijKL, KLMN >( einsum< ijKL, IJKL >( dTau_dPK2, dPk2_dC ), dC_dF ) + dTau_dF;
  }

  StateView ViscoelasticNeoHooke::getStateView( const std::string& stateName )
  {
    return stateVars->getStateView( stateName );
  }

  void ViscoelasticNeoHooke::assignStateVars( double* stateVars_, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                                << ": Not sufficient "
                                                   "stateVars!" );

    this->stateVars = std::make_unique< ViscoelasticNeoHookeStateVarManager >( stateVars_ );
  }
} // namespace Marmot::Materials
