#include "Marmot/MarmotFiniteStrainViscoelasticity.h"

namespace Marmot {
  namespace ContinuumMechanics::FiniteStrain::Viscoelasticity {

    MaxwellProperties createMaxwellProperties( int nMaxwell, const double* gammaTauPairVector )
    {

      std::vector< double > gamma( nMaxwell );
      std::vector< double > tau( nMaxwell );
      double                sumGamma = 0.0;
      for ( int i = 0; i < nMaxwell; ++i ) {
        gamma[i] = gammaTauPairVector[i * 2];
        tau[i]   = gammaTauPairVector[i * 2 + 1];
        sumGamma += gamma[i];
      }
      return MaxwellProperties( nMaxwell, gamma, sumGamma, tau );
    }

    void evaluateGeneralizedMaxwellModel( FastorStandardTensors::Tensor33d&           stress,
                                          FastorStandardTensors::Tensor3333d&         tangent,
                                          const FastorStandardTensors::Tensor333333d& dTangent_dDeformation,
                                          const FastorStandardTensors::Tensor3333d&   initialCompliance,
                                          const FastorStandardTensors::Tensor33d&     dStress,
                                          const double                                dT,
                                          const MaxwellProperties&                    maxwellProperties,
                                          double*                                     stateVars )
    {

      if ( maxwellProperties.nMaxwell == 0 )
        return;

      using namespace Fastor;
      using namespace Marmot::FastorStandardTensors;
      using namespace Marmot::FastorIndices;

      // copy of tangent to be incremented
      Tensor3333d initialTangent = tangent;

      // scale equilibrium stress contribution
      stress *= ( 1.0 - maxwellProperties.sumGamma );
      tangent *= ( 1.0 - maxwellProperties.sumGamma );

      for ( size_t i = 0; i < maxwellProperties.nMaxwell; ++i ) {
        // get old  maxewell element stress from state variables
        const Tensor33d& H_n = Tensor33d( stateVars + i * 9 );

        // get parameters of maxwell element
        const double& tau   = maxwellProperties.tau[i];
        const double& gamma = maxwellProperties.gamma[i];

        const double dT_tau    = std::max( dT / tau, 1e-15 );
        const double expFactor = Math::exp( -dT_tau );

        // compute new stress in maxwell element
        const Tensor33d H_np = expFactor * H_n + gamma / dT_tau * ( 1.0 - expFactor ) * dStress;

        // add contribution to stress
        const Tensor33d dStress_ = einsum< ijkl, ij >( initialTangent, einsum< ijkl, ij >( initialCompliance, H_np ) );
        stress += dStress_;

        // update tangent
        tangent += gamma / dT_tau * ( 1.0 - expFactor ) *
                   einsum< ijkl, mnij >( initialTangent, einsum< ijkl, mnij >( initialTangent, initialCompliance ) );
        tangent += einsum< ijklmn, ij >( dTangent_dDeformation, einsum< ijkl, kl >( initialCompliance, H_np ) );

        // update state variables
        memcpy( stateVars + i * 9, &H_np, 9 * sizeof( double ) );
      }
    }

  } // namespace ContinuumMechanics::FiniteStrain::Viscoelasticity
} // namespace Marmot
