#include "Marmot/FiniteStrainIsotropicBiotViscoelasticity.h"
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEigenSystems.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
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

  FiniteStrainIsotropicBiotViscoelasticity::FiniteStrainIsotropicBiotViscoelasticity( const double* materialProperties,
                                                                                      int           nMaterialProperties,
                                                                                      int           materialLabel )
    : MarmotMaterialFiniteStrainAD( materialProperties, nMaterialProperties, materialLabel ),
      K( materialProperties[0] ),
      G( materialProperties[1] ),
      maxwellProperties(
        ContinuumMechanics::FiniteStrain::Viscoelasticity::createMaxwellProperties( materialProperties[2],
                                                                                    &materialProperties[3] ) )
  {

    Tensor3333d initialStiffness = std::get< 2 >(
      ContinuumMechanics::EnergyDensityFunctions::SecondOrderDerived::BiotNeoHooke( evaluate( Spatial3D::I ), K, G ) );

    Fastor::Tensor< double, 6, 6 > reducedStiffness;
    Fastor::Tensor< double, 9, 9 > initialStiffness99 = Fastor::reshape< 9, 9 >( initialStiffness );
    initialStiffness99 = 0.5 * ( initialStiffness99 + transpose( initialStiffness99 ) ); // make sure it's symmetric
    int rowColMap[6]   = { 0, 1, 2, 4, 5, 8 };
    // double                         scaling[6]         = { 1.0, 2.0, 2.0, 1.0, 2.0, 1.0 };
    int rowColMapFull[9] = { 0, 1, 2, 1, 3, 4, 2, 4, 5 };

    // manually modify
    initialStiffness99( 1, all ) += initialStiffness99( 3, all );
    initialStiffness99( 3, all ) = initialStiffness99( 1, all );
    initialStiffness99( 2, all ) += initialStiffness99( 6, all );
    initialStiffness99( 6, all ) = initialStiffness99( 2, all );
    initialStiffness99( 5, all ) += initialStiffness99( 7, all );
    initialStiffness99( 7, all ) = initialStiffness99( 5, all );

    for ( int i = 0; i < 6; ++i )
      for ( int j = 0; j < 6; ++j )
        reducedStiffness( i, j ) = initialStiffness99( rowColMap[i], rowColMap[j] );

    // invert reduced stiffness
    Fastor::Tensor< double, 6, 6 > reducedCompliance = inverse< InvCompType::BlockLU >( reducedStiffness );

    Fastor::Tensor< double, 9, 9 > fullCompliance( 0.0 );
    // size_t                         rowColMapFull[9] = { 0, 1, 2, 1, 3, 4, 2, 4, 5 };
    double scalingFull[9] = { 1, 0.5, 0.5, .5, 1, .5, .5, .5, 1 };

    for ( int i = 0; i < 9; ++i )
      for ( int j = 0; j < 9; ++j )
        fullCompliance( i, j ) = reducedCompliance( rowColMapFull[i], rowColMapFull[j] ) * scalingFull[j];

    // reshape back to 4th order
    Tensor3333d initialComplianceReal = Fastor::reshape< 3, 3, 3, 3 >( fullCompliance );

    if ( norm( initialComplianceReal ) != norm( initialComplianceReal ) ) {
      std::cout << "Initial Compliance: \n" << initialComplianceReal << std::endl;
      std::cout << "Initial Stiffness: \n" << initialStiffness99 << std::endl;
      throw std::runtime_error(
        "CompressibleFiniteStrainLinearViscoelasticity::CompressibleFiniteStrainLinearViscoelasticity: Initial "
        "compliance has NaN values. Check material properties." );
    }
    initialCompliance = makeDual( initialComplianceReal );

    initializeStateLayout();
  }

  void FiniteStrainIsotropicBiotViscoelasticity::computeStressAD( ConstitutiveResponseAD< 3 >& response,
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

    auto [psi, S_biot, dS_biot_dU] = EnergyDensityFunctions::SecondOrderDerived::BiotNeoHooke( U, K, G );

    Tensor33d& S_biot_old = stateLayout.getAs< Tensor33d& >( response.stateVars, "S0_old" );

    const Tensor33t< scalar > dS_biot = S_biot - makeDual( S_biot_old );
    memcpy( S_biot_old.data(), makeReal( S_biot ).data(), 9 * sizeof( double ) );

    // add viscoelastic contribution to deviatoric PK2 stress
    ContinuumMechanics::FiniteStrain::Viscoelasticity::evaluateGeneralizedMaxwellModel<
      scalar >( S_biot,
                dS_biot_dU,
                initialCompliance,
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
    response.elasticEnergyDensity = Math::makeReal( psi );
  }
} // namespace Marmot::Materials
