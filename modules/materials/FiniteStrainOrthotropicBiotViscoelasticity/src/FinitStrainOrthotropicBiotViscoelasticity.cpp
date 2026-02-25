#include "Marmot/FiniteStrainOrthotropicBiotViscoelasticity.h"
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEigenSystems.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteStrainViscoelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
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

  FiniteStrainOrthotropicBiotViscoelasticity::FiniteStrainOrthotropicBiotViscoelasticity(
    const double* materialProperties,
    int           nMaterialProperties,
    int           materialLabel )
    : MarmotMaterialFiniteStrainAD( materialProperties, nMaterialProperties, materialLabel ),
      E1( materialProperties[0] ),
      E2( materialProperties[1] ),
      E3( materialProperties[2] ),
      nu12( materialProperties[3] ),
      nu13( materialProperties[4] ),
      nu23( materialProperties[5] ),
      G12( materialProperties[6] ),
      G13( materialProperties[7] ),
      G23( materialProperties[8] ),
      maxwellProperties(
        ContinuumMechanics::FiniteStrain::Viscoelasticity::createMaxwellProperties( materialProperties[9],
                                                                                    &materialProperties[10] ) )
  {

    Matrix6d C_voigt = ContinuumMechanics::Elasticity::Orthotropic::stiffnessTensor( E1,
                                                                                     E2,
                                                                                     E3,
                                                                                     nu12,
                                                                                     nu23,
                                                                                     nu13,
                                                                                     G12,
                                                                                     G23,
                                                                                     G13 );

    // populate dBiotStress_dU tensor from voigt matrix
    // Voigt notation: 11->0, 22->1, 33->2, 23->3, 13->4, 12->5
    dBiotStress_dU( 0, 0, 0, 0 ) = C_voigt( 0, 0 );
    dBiotStress_dU( 1, 1, 1, 1 ) = C_voigt( 1, 1 );
    dBiotStress_dU( 2, 2, 2, 2 ) = C_voigt( 2, 2 );
    dBiotStress_dU( 1, 1, 0, 0 ) = dBiotStress_dU( 0, 0, 1, 1 ) = C_voigt( 1, 0 );
    dBiotStress_dU( 2, 2, 0, 0 ) = dBiotStress_dU( 0, 0, 2, 2 ) = C_voigt( 2, 0 );
    dBiotStress_dU( 2, 2, 1, 1 ) = dBiotStress_dU( 1, 1, 2, 2 ) = C_voigt( 2, 1 );
    dBiotStress_dU( 0, 1, 0, 1 ) = dBiotStress_dU( 1,
                                                   0,
                                                   0,
                                                   1 ) = dBiotStress_dU( 0,
                                                                         1,
                                                                         1,
                                                                         0 ) = dBiotStress_dU( 1,
                                                                                               0,
                                                                                               1,
                                                                                               0 ) = C_voigt( 3, 3 );
    dBiotStress_dU( 1, 2, 1, 2 ) = dBiotStress_dU( 2,
                                                   1,
                                                   1,
                                                   2 ) = dBiotStress_dU( 1,
                                                                         2,
                                                                         2,
                                                                         1 ) = dBiotStress_dU( 2,
                                                                                               1,
                                                                                               2,
                                                                                               1 ) = C_voigt( 5, 5 );
    dBiotStress_dU( 0, 2, 0, 2 ) = dBiotStress_dU( 2,
                                                   0,
                                                   0,
                                                   2 ) = dBiotStress_dU( 0,
                                                                         2,
                                                                         2,
                                                                         0 ) = dBiotStress_dU( 2,
                                                                                               0,
                                                                                               2,
                                                                                               0 ) = C_voigt( 4, 4 );

    initializeStateLayout();
  }

  void FiniteStrainOrthotropicBiotViscoelasticity::computeStressAD( ConstitutiveResponseAD< 3 >& response,
                                                                    const DeformationAD< 3 >&    deformation,
                                                                    const TimeIncrement&         timeIncrement ) const
  {

    using scalar  = autodiff::dual;
    const auto& F = deformation.F;

    using namespace ContinuumMechanics;
    // compute Cauchy-Green deformation
    const Tensor33t< scalar > C = DeformationMeasures::rightCauchyGreen( F );

    auto [lam, Q] = Math::computeEigenSystemJacobi( C );
    Tensor33t< scalar > principalStretch( 0. );
    for ( int i = 0; i < 3; ++i ) {
      principalStretch( i, i ) = sqrt( lam( i ) );
    }
    const Tensor33t< scalar > U = Q % principalStretch % transpose( Q );

    const Tensor33t< scalar > I = makeDual( Spatial3D::I );

    Tensor33t< scalar > S_biot = einsum< ijkl, kl >( makeDual( dBiotStress_dU ), evaluate( U - I ) );

    Tensor33d& S_biot_old = stateLayout.getAs< Tensor33d& >( response.stateVars, "S0_old" );

    const Tensor33t< scalar > dS_biot = S_biot - makeDual( S_biot_old );
    memcpy( S_biot_old.data(), makeReal( S_biot ).data(), 9 * sizeof( double ) );

    // add viscoelastic contribution to deviatoric PK2 stress
    ContinuumMechanics::FiniteStrain::Viscoelasticity::evaluateGeneralizedMaxwellModel<
      scalar >( S_biot,
                dS_biot,
                timeIncrement.dT,
                maxwellProperties,
                stateLayout.getPtr( response.stateVars, "creepStateVars" ) );

    Tensor33t< scalar > S_biot_rotated = transpose( Q ) % S_biot % Q;
    Tensor33t< scalar > PK2_rotated( 0 );

    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        PK2_rotated( i, j ) = 2.0 * S_biot_rotated( i, j ) / ( principalStretch( i, i ) + principalStretch( j, j ) );
      }
    }
    Tensor33t< scalar > PK2 = Q % PK2_rotated % transpose( Q );

    const Tensor33t< scalar > tau = StressMeasures::KirchhoffStressFromPK2( PK2, F );
    response.tau                  = tau;
    response.rho                  = 1.0;
    response.elasticEnergyDensity = -1;
  }
} // namespace Marmot::Materials
