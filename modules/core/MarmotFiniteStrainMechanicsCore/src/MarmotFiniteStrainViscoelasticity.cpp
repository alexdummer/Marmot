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
      const Tensor3333d initialTangent = tangent;

      // scale equilibrium stress contribution
      stress  = stress * ( 1.0 - maxwellProperties.sumGamma );
      tangent = tangent * ( 1.0 - maxwellProperties.sumGamma );

      for ( size_t i = 0; i < maxwellProperties.nMaxwell; ++i ) {
        // get old  maxewell element stress from state variables
        const Tensor33d& Q_n = Tensor33d( stateVars + i * 9 );

        // get parameters of maxwell element
        const double& tau   = maxwellProperties.tau[i];
        const double& gamma = maxwellProperties.gamma[i];

        const double dT_tau    = std::max( dT / tau, 1e-15 );
        const double expFactor = Math::exp( -dT_tau );

        double alpha = expFactor;
        double beta  = gamma / dT_tau * ( 1.0 - expFactor );

        if ( dT_tau < 1e-6 ) {
          // use taylor expansion for small dt/tau
          alpha = 1.0 - dT_tau + 0.5 * dT_tau * dT_tau;
          beta  = gamma * ( 1.0 - 0.5 * dT_tau + 1.0 / 6.0 * dT_tau * dT_tau );
        }

        // compute new stress in maxwell element
        const Tensor33d Q_np = alpha * Q_n + beta * dStress;

        // add contribution to stress
        const Tensor33d H_np = einsum< ijkl, kl >( initialCompliance, Q_np );
        stress += einsum< ij, ijkl >( H_np, initialTangent );

        const Tensor3333d dH_np_dDeformation = einsum< ijmn, mnKL >( initialCompliance,
                                                                     evaluate( beta * initialTangent ) );

        tangent += einsum< ijkl, ijmn >( initialTangent, dH_np_dDeformation );
        tangent += einsum< ij, ijklmn >( H_np, dTangent_dDeformation );

        // update state variables
        memcpy( stateVars + i * 9, &Q_np, 9 * sizeof( double ) );
      }
    }

  } // namespace ContinuumMechanics::FiniteStrain::Viscoelasticity
} // namespace Marmot
