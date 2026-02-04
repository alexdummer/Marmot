#include "Marmot/CompressibleFiniteStrainLinearViscoelasticity.h"
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotDeformationMeasures.h"
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

  std::tuple< double,
              FastorStandardTensors::Tensor33d,
              FastorStandardTensors::Tensor3333d,
              FastorStandardTensors::Tensor333333d >
  CompressibleFiniteStrainLinearViscoelasticity::computeEnergyDensityAndDerivatives(
    const FastorStandardTensors::Tensor33d& C ) const
  {

    // select hyperelastic base
    // compute enregy density and derivatives
    // with autodiff
    // we need a function mapping from tensor33 with dual3rd  to dual3rd
    std::function< autodiff::dual3rd( const FastorStandardTensors::Tensor33t< autodiff::dual3rd >& ) >
      energyDensityFunction;

    switch ( hyperelasticBase ) {
    case NeoHooke:
      return ContinuumMechanics::EnergyDensityFunctions::ThirdOrderDerived::standardNeoHooke( C,
                                                                                              elasticProperties[0],
                                                                                              elasticProperties[1] );
      break;
    case Yeoh:
      energyDensityFunction = [this]( const FastorStandardTensors::Tensor33t< autodiff::dual3rd >& C_ad ) {
        autodiff::dual3rd
          res = ContinuumMechanics::EnergyDensityFunctions::YeohPotential< autodiff::dual3rd >( C_ad,
                                                                                                elasticProperties[0],
                                                                                                elasticProperties[1],
                                                                                                elasticProperties[2],
                                                                                                elasticProperties[3] );
        return res;
      };
      break;
    case MooneyRivlin:
      energyDensityFunction = [this]( const FastorStandardTensors::Tensor33t< autodiff::dual3rd >& C_ad ) {
        autodiff::dual3rd res = ContinuumMechanics::EnergyDensityFunctions::MooneyRivlinPotential<
          autodiff::dual3rd >( C_ad, elasticProperties[0], elasticProperties[1], elasticProperties[2] );
        return res;
      };
      break;
    default:
      throw std::runtime_error( "CompressibleFiniteStrainLinearViscoelasticity::computeEnergyDensityAndDerivatives: "
                                "Unknown hyperelastic base." );
    }
    return Marmot::AutomaticDifferentiation::ThirdOrder::d3f_dT3< 3 >( energyDensityFunction, C );
  }

  CompressibleFiniteStrainLinearViscoelasticity::CompressibleFiniteStrainLinearViscoelasticity(
    const double* materialProperties,
    int           nMaterialProperties,
    int           materialLabel )
    : MarmotMaterialFiniteStrain( materialProperties, nMaterialProperties, materialLabel ),
      hyperelasticBase( static_cast< HyperelasticBase >( materialProperties[0] ) ),
      onlyShearCreep( materialProperties[1] ),
      elasticProperties( &materialProperties[2], nElasticPropertiesMap.at( hyperelasticBase ) ),
      maxwellProperties(
        ContinuumMechanics::FiniteStrain::Viscoelasticity::
          createMaxwellProperties( materialProperties[2 + nElasticPropertiesMap.at( hyperelasticBase )],
                                   &materialProperties[3 + nElasticPropertiesMap.at( hyperelasticBase )] ) )
  {

    // compute initial deviatoric compliance tensor
    // start from initial stiffness tensor
    Tensor3333d initialStiffness = 4. * std::get< 2 >( computeEnergyDensityAndDerivatives( evaluate( Spatial3D::I ) ) );

    // initialStiffness = initialStiffness - initialStiffness % Spatial3D::I;
    // initialStiffness = initialStiffness - ( 1 / 3.0 ) * outer( Spatial3D::I, einsum< ijkl, kl >( initialStiffness,
    // Spatial3D::I ) );

    Fastor::Tensor< double, 6, 6 > reducedStiffness;
    Fastor::Tensor< double, 9, 9 > initialStiffness99 = Fastor::reshape< 9, 9 >( initialStiffness );
    initialStiffness99  = 0.5 * ( initialStiffness99 + transpose( initialStiffness99 ) ); // make sure it's symmetric
    size_t rowColMap[6] = { 0, 1, 2, 4, 5, 8 };
    // double                         scaling[6]         = { 1.0, 2.0, 2.0, 1.0, 2.0, 1.0 };
    size_t rowColMapFull[9] = { 0, 1, 2, 1, 3, 4, 2, 4, 5 };

    // manually modify
    initialStiffness99( 1, all ) += initialStiffness99( 3, all );
    initialStiffness99( 3, all ) = initialStiffness99( 1, all );
    initialStiffness99( 2, all ) += initialStiffness99( 6, all );
    initialStiffness99( 6, all ) = initialStiffness99( 2, all );
    initialStiffness99( 5, all ) += initialStiffness99( 7, all );
    initialStiffness99( 7, all ) = initialStiffness99( 5, all );

    for ( size_t i = 0; i < 6; ++i )
      for ( size_t j = 0; j < 6; ++j )
        reducedStiffness( i, j ) = initialStiffness99( rowColMap[i], rowColMap[j] );

    // invert reduced stiffness
    Fastor::Tensor< double, 6, 6 > reducedCompliance = inverse< InvCompType::BlockLU >( reducedStiffness );

    Fastor::Tensor< double, 9, 9 > fullCompliance( 0.0 );
    // size_t                         rowColMapFull[9] = { 0, 1, 2, 1, 3, 4, 2, 4, 5 };
    double scalingFull[9] = { 1, 0.5, 0.5, .5, 1, .5, .5, .5, 1 };

    for ( size_t i = 0; i < 9; ++i )
      for ( size_t j = 0; j < 9; ++j )
        fullCompliance( i, j ) = reducedCompliance( rowColMapFull[i], rowColMapFull[j] ) * scalingFull[j];

    // reshape back to 4th order
    initialCompliance = Fastor::reshape< 3, 3, 3, 3 >( fullCompliance );

    if ( norm( initialCompliance ) != norm( initialCompliance ) ) {
      std::cout << "Initial Compliance: \n" << initialCompliance << std::endl;
      std::cout << "Initial Stiffness: \n" << initialStiffness99 << std::endl;
      throw std::runtime_error(
        "CompressibleFiniteStrainLinearViscoelasticity::CompressibleFiniteStrainLinearViscoelasticity: Initial "
        "compliance has NaN values. Check material properties." );
    }

    initializeStateLayout();
  }

  void CompressibleFiniteStrainLinearViscoelasticity::computeStress( ConstitutiveResponse< 3 >& response,
                                                                     AlgorithmicModuli< 3 >&    tangents,
                                                                     const Deformation< 3 >&    deformation,
                                                                     const TimeIncrement&       timeIncrement ) const
  {
    const auto& F = deformation.F;

    using namespace ContinuumMechanics;
    // compute Cauchy-Green deformation
    const auto [C, dC_dF] = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( F );

    // compute energy density, first, second and third partial derivatives wrt Cauchy Green deformation
    const auto [psi_, dPsi_dC, d2Psi_dCdC, d3Psi_dCdCdC] = computeEnergyDensityAndDerivatives( C );

    // compute initial PK2 stress
    const Tensor33d PK2zero = 2. * dPsi_dC;

    // split into volumetric and deviatoric parts
    // we only use the deviatoric part for viscoelastic evolution
    const Tensor33d PK2vol = ( onlyShearCreep / 3.0 ) * trace( PK2zero ) * Spatial3D::I;
    Tensor33d       PK2dev = PK2zero - PK2vol;

    // conpute increment in initial PK2dev
    Tensor33d&      PK2dev_old = stateLayout.getAs< Tensor33d& >( response.stateVars, "S0_old" );
    const Tensor33d dPK2       = PK2dev - PK2dev_old;

    const Tensor3333d dPK2_dE    = 4. * d2Psi_dCdC;
    const Tensor3333d dPK2vol_dE = onlyShearCreep / 3.0 *
                                   outer( Spatial3D::I, einsum< ijkl, kl >( dPK2_dE, Spatial3D::I ) );
    Tensor3333d         dPK2dev_dE    = 4. * d2Psi_dCdC - dPK2vol_dE;
    const Tensor333333d d2PK2dev_dEdE = 8. * d3Psi_dCdCdC;

    // add viscoelastic contribution to deviatoric PK2 stress
    ContinuumMechanics::FiniteStrain::Viscoelasticity::evaluateGeneralizedMaxwellModel( PK2dev,
                                                                                        dPK2dev_dE,
                                                                                        d2PK2dev_dEdE,
                                                                                        initialCompliance,
                                                                                        dPK2,
                                                                                        timeIncrement.dT,
                                                                                        maxwellProperties,
                                                                                        stateLayout
                                                                                          .getPtr( response.stateVars,
                                                                                                   "creepStateVars" ) );

    const Tensor33d PK2 = evaluate( PK2vol + PK2dev );

    // compute Kirchhoff stress
    const auto [tau, dTau_dPK2, dTau_dF] = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, F );
    response.tau                         = tau;
    response.rho                         = 1.0;
    response.elasticEnergyDensity        = psi_;

    // compute tangent operator
    tangents.dTau_dF = einsum< ijKL, KLMN >( einsum< ijKL, KLMN >( dTau_dPK2, evaluate( dPK2vol_dE + dPK2dev_dE ) ),
                                             0.5 * dC_dF ) +
                       dTau_dF;

    memcpy( PK2dev_old.data(), PK2dev.data(), 9 * sizeof( double ) );
  }
} // namespace Marmot::Materials
