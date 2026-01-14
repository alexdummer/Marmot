#include "Marmot/LinearViscoelasticCompressibleMooneyRivlin.h"
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

  // use autodiff to get stress and tangent
  std::tuple< double, Tensor33d, Tensor3333d, Tensor333333d > computeStressAndTangentAutodiff( const Tensor33d& C,
                                                                                               const double     C1,
                                                                                               const double     C2,
                                                                                               const double     K )
  {

    auto func = [&]( const Tensor33t< autodiff::dual3rd >& C_var ) {
      autodiff::dual3rd res = ContinuumMechanics::EnergyDensityFunctions::MooneyRivlinPotential( C_var, C1, C2, 2 / K );
      return res;
    };

    auto [psi, dPsi_dC, d2Psi_dCdC, d3Psi_dCdCdC] = Marmot::AutomaticDifferentiation::ThirdOrder::d3f_dT3< 3 >( func,
                                                                                                                C );

    return { psi, dPsi_dC, d2Psi_dCdC, d3Psi_dCdCdC };
  }
  LinearViscoelasticCompressibleMooneyRivlin::LinearViscoelasticCompressibleMooneyRivlin(
    const double* materialProperties,
    int           nMaterialProperties,
    int           materialLabel )
    : MarmotMaterialFiniteStrain( materialProperties, nMaterialProperties, materialLabel ),
      C1( materialProperties[0] ),
      C2( materialProperties[1] ),
      K( materialProperties[2] ),
      maxwellProperties(
        ContinuumMechanics::FiniteStrain::Viscoelasticity::createMaxwellProperties( materialProperties[3],
                                                                                    &materialProperties[4] ) )
  {

    // compute initial compliance
    auto I = Spatial3D::I;

    // print parameters
    std::cout << "Material Parameters: " << std::endl;
    std::cout << "C1: " << C1 << std::endl;
    std::cout << "C2: " << C2 << std::endl;
    std::cout << "K: " << K << std::endl;

    std::cout << "D1: " << 2. / K << std::endl;

    auto [psi_init, dPsi_dC_init, d2Psi_dCdC_init, d3Psi_dCdCdC_init] = computeStressAndTangentAutodiff( I, C1, C2, K );

    // print the initial stiffness
    std::cout << "Initial Energy Density: " << psi_init << std::endl;
    std::cout << "Initial First Derivative dPsi/dC: \n" << dPsi_dC_init << std::endl;
    // std::cout << "Initial Second Derivative d2Psi/dCdC: \n" << d2Psi_dCdC_init << std::endl;
    // std::cout << "Initial Third Derivative d3Psi/dCdCdC: \n" << d3Psi_dCdCdC_init << std::endl;

    Tensor3333d initialStiffness = 4 * d2Psi_dCdC_init;
    // - ( 4.0 / 3.0 ) * outer( Spatial3D::I, einsum< ijkl, kl >( d2Psi_dCdC_init, Spatial3D::I ) ) ;

    // Fastor::Tensor< double, 3, 3, 3, 3 > initialStiffness = 2. *
    //     std::get<2>( ContinuumMechanics::EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB(Spatial3D::I,
    //     K, G ) );

    std::cout << "Initial Stiffness Tensor: \n" << initialStiffness << std::endl;

    auto Isymm = Spatial3D::ISymm;

    Fastor::Tensor< double, 6, 6 > reducedStiffness;
    Fastor::Tensor< double, 9, 9 > initialStiffness99 = Fastor::reshape< 9, 9 >(
      einsum< ijkl, ijmn >( initialStiffness, Spatial3D::ISymm ) );
    initialStiffness99 = 0.5 * ( initialStiffness99 + transpose( initialStiffness99 ) ); // make sure it's symmetric
    std::cout << "Initial Stiffness 9x9: \n" << initialStiffness99 << std::endl;
    size_t rowColMap[6] = { 0, 1, 2, 4, 5, 8 };
    double scaling[6]   = { 1.0, 2.0, 2.0, 1.0, 2.0, 1.0 };
    for ( size_t i = 0; i < 6; ++i )
      for ( size_t j = 0; j < 6; ++j )
        reducedStiffness( i, j ) = initialStiffness99( rowColMap[i], rowColMap[j] ) * scaling[j];

    std::cout << "Reduced Stiffness Tensor: \n" << reducedStiffness << std::endl;

    // invert reduced stiffness
    Fastor::Tensor< double, 6, 6 > reducedCompliance = inverse< InvCompType::BlockLU >( reducedStiffness );

    // check if inversion is correct
    Fastor::Tensor< double, 6, 6 > checkIdentity6 = reducedStiffness % reducedCompliance;
    std::cout << "Check Identity (6x6): \n" << checkIdentity6 << std::endl;

    Fastor::Tensor< double, 9, 9 > fullCompliance( 0.0 );
    size_t                         rowColMapFull[9] = { 0, 1, 2, 1, 3, 4, 2, 4, 5 };
    double                         scalingFull[9]   = { 1, 0.5, 0.5, .5, 1, .5, .5, .5, 1 };

    for ( size_t i = 0; i < 9; ++i )
      for ( size_t j = 0; j < 9; ++j )
        fullCompliance( i, j ) = reducedCompliance( rowColMapFull[i], rowColMapFull[j] ) * scalingFull[j];

    // reshape back to 4th order
    initialCompliance = Fastor::reshape< 3, 3, 3, 3 >( fullCompliance );

    // std::cout << "Initial Compliance Tensor: \n" << initialCompliance << std::endl;
    // std::exit(0);

    // check if inversion is correct
    Fastor::Tensor< double, 9, 9 > checkIdentity = initialStiffness99 % fullCompliance;
    std::cout << "Check Identity (9x9): \n" << checkIdentity << std::endl;

    initializeStateLayout();

    // std::exit(0);
  }

  void LinearViscoelasticCompressibleMooneyRivlin::computeStress( ConstitutiveResponse< 3 >& response,
                                                                  AlgorithmicModuli< 3 >&    tangents,
                                                                  const Deformation< 3 >&    deformation,
                                                                  const TimeIncrement&       timeIncrement ) const
  {
    const auto& F_ = deformation.F;

    using namespace ContinuumMechanics;
    // compute Cauchy-Green deformation
    const auto [C, dC_dF] = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( F_ );

    // compute energy density, first and second partial derivatives wrt Cauchy Green deformation
    const auto [psi_, dPsi_dC, d2Psi_dCdC, d3Psi_dCdCdC] = computeStressAndTangentAutodiff( C, C1, C2, K );

    // compute Kirchhoff stress
    Tensor33d PK2 = 2. * dPsi_dC;

    // split into volumetric and deviatoric parts
    Tensor33d PK2vol = ( 1.0 / 3.0 ) * trace( PK2 ) * Spatial3D::I;
    Tensor33d PK2dev = PK2 - PK2vol;

    // delta PK2dev
    Tensor33d dPK2 = PK2dev - stateLayout.getAs< Tensor33d& >( response.stateVars, "S0_old" );
    memcpy( stateLayout.getPtr( response.stateVars, "S0_old" ), PK2dev.data(), 9 * sizeof( double ) );

    Tensor3333d   dPK2_dE       = 4. * d2Psi_dCdC;
    Tensor3333d   dPK2vol_dE    = 1. / 3.0 * outer( Spatial3D::I, einsum< ijkl, kl >( dPK2_dE, Spatial3D::I ) );
    Tensor3333d   dPK2dev_dE    = 4. * d2Psi_dCdC - dPK2vol_dE;
    Tensor333333d d2PK2dev_dEdE = 8. * d3Psi_dCdCdC;

    // add viscoelastic contribution to PK2 stress
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
    // compute Kirchhoff stress
    PK2 = PK2vol + PK2dev;

    const auto [tau, dTau_dPK2, dTau_dF] = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, F_ );
    response.tau                         = tau;
    response.rho                         = 1.0;
    response.elasticEnergyDensity        = psi_;

    // compute tangent operator
    tangents.dTau_dF = einsum< ijKL, KLMN >( einsum< ijKL, IJKL >( dTau_dPK2, dPK2vol_dE + dPK2dev_dE ), 0.5 * dC_dF ) +
                       dTau_dF;
  }
} // namespace Marmot::Materials
