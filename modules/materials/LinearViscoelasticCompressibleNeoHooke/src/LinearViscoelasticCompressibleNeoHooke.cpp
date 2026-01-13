#include "Marmot/LinearViscoelasticCompressibleNeoHooke.h"
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

  LinearViscoelasticCompressibleNeoHooke::LinearViscoelasticCompressibleNeoHooke( const double* materialProperties,
                                                                                  int           nMaterialProperties,
                                                                                  int           materialLabel )
    : MarmotMaterialFiniteStrain( materialProperties, nMaterialProperties, materialLabel ),
      K( materialProperties[0] ),
      G( materialProperties[1] ),
      maxwellProperties(
        ContinuumMechanics::FiniteStrain::Viscoelasticity::createMaxwellProperties( materialProperties[2],
                                                                                    &materialProperties[3] ) )
  {

    // compute initial compliance

    Tensor3333d initialStiffness = 4. * std::get< 2 >( ContinuumMechanics::EnergyDensityFunctions::ThirdOrderDerived::
                                                         standardNeoHooke( Spatial3D::I, K, G ) );

    // Fastor::Tensor< double, 3, 3, 3, 3 > initialStiffness = 2. *
    //     std::get<2>( ContinuumMechanics::EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB(Spatial3D::I,
    //     K, G ) );

    std::cout << "Initial Stiffness Tensor: \n" << initialStiffness << std::endl;

    Fastor::Tensor< double, 6, 6 > reducedStiffness;
    Fastor::Tensor< double, 9, 9 > initialStiffness99 = Fastor::reshape< 9, 9 >( initialStiffness );
    std::cout << "Initial Stiffness 9x9: \n" << initialStiffness99 << std::endl;
    size_t rowColMap[6] = { 0, 1, 2, 4, 5, 8 };
    double scaling[6]   = { 1.0, 2.0, 2.0, 1.0, 2.0, 1.0 };
    for ( size_t i = 0; i < 6; ++i )
      for ( size_t j = 0; j < 6; ++j )
        reducedStiffness( i, j ) = initialStiffness99( rowColMap[i], rowColMap[j] ) * scaling[j];

    std::cout << "Reduced Stiffness Tensor: \n" << reducedStiffness << std::endl;

    // invert reduced stiffness
    Fastor::Tensor< double, 6, 6 > reducedCompliance = inverse< InvCompType::BlockLU >( reducedStiffness );

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

    initializeStateLayout();

    // std::exit(0);
  }

  void LinearViscoelasticCompressibleNeoHooke::computeStress( ConstitutiveResponse< 3 >& response,
                                                              AlgorithmicModuli< 3 >&    tangents,
                                                              const Deformation< 3 >&    deformation,
                                                              const TimeIncrement&       timeIncrement ) const
  {
    const auto& F_ = deformation.F;

    using namespace ContinuumMechanics;
    // compute Cauchy-Green deformation
    const auto [C, dC_dF] = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( F_ );

    // compute energy density, first and second partial derivatives wrt Cauchy Green deformation
    const auto [psi_,
                dPsi_dC,
                d2Psi_dCdC,
                d3Psi_dCdCdC] = EnergyDensityFunctions::ThirdOrderDerived::standardNeoHooke( C, K, G );

    // compute Kirchhoff stress
    Tensor33d PK2  = 2. * dPsi_dC;
    Tensor33d dPK2 = PK2 - stateLayout.getAs< Tensor33d& >( response.stateVars, "S0_old" );
    memcpy( stateLayout.getPtr( response.stateVars, "S0_old" ), PK2.data(), 9 * sizeof( double ) );

    Tensor3333d   dPK2_dE    = 4. * d2Psi_dCdC;
    Tensor333333d d2PK2_dEdE = 8. * d3Psi_dCdCdC;

    // add viscoelastic contribution to PK2 stress
    ContinuumMechanics::FiniteStrain::Viscoelasticity::evaluateGeneralizedMaxwellModel( PK2,
                                                                                        dPK2_dE,
                                                                                        d2PK2_dEdE,
                                                                                        initialCompliance,
                                                                                        dPK2,
                                                                                        timeIncrement.dT,
                                                                                        maxwellProperties,
                                                                                        stateLayout
                                                                                          .getPtr( response.stateVars,
                                                                                                   "creepStateVars" ) );

    const auto [tau, dTau_dPK2, dTau_dF] = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, F_ );
    response.tau                         = tau;
    response.rho                         = 1.0;
    response.elasticEnergyDensity        = psi_;

    // compute tangent operator
    tangents.dTau_dF = einsum< ijKL, KLMN >( einsum< ijKL, IJKL >( dTau_dPK2, dPK2_dE ), 0.5 * dC_dF ) + dTau_dF;
  }
} // namespace Marmot::Materials
