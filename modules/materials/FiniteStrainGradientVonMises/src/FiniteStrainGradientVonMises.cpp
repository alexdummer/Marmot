/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

#include "Marmot/FiniteStrainGradientVonMises.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteStrainPlasticity.h"
#include "Marmot/MarmotStressMeasures.h"

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace FastorIndices;
  using namespace FastorStandardTensors;
  using namespace autodiff;
  using namespace ContinuumMechanics;

  FiniteStrainGradientVonMises::FiniteStrainGradientVonMises( const double* materialProperties,
                                                              int           nMaterialProperties,
                                                              int           materialNumber )
    : MarmotMaterialGradientPlasticityFiniteStrainAD< 1 >( materialProperties, nMaterialProperties, materialNumber ),
      K( materialProperties[0] ),
      G( materialProperties[1] ),
      fy0( materialProperties[2] ),
      H( materialProperties[3] ),
      g( materialProperties[4] )
  {
    stateLayout.add( "Fp", 9 );
    stateLayout.add( "kappa", 1 );
    stateLayout.add( "laplaceKappa", 1 );
    stateLayout.finalize();
  }

  double FiniteStrainGradientVonMises::getDensity( const double* stateVars ) const
  {
    if ( nMaterialProperties <= 5 )
      throw std::runtime_error( "Density not provided in material properties." );
    return materialProperties[5];
  }

  std::vector< double > FiniteStrainGradientVonMises::getNonlocalViscosity( const double* stateVars ) const
  {
    if ( nMaterialProperties <= 6 )
      throw std::runtime_error( "Viscosity not provided in material properties." );
    return { materialProperties[6] };
  }

  void FiniteStrainGradientVonMises::initializeYourself( double* stateVars, int nStateVars )
  {
    for ( int i = 0; i < nStateVars; ++i ) {
      stateVars[i] = 0.0;
    }

    TensorMap33d Fp = stateLayout.getAs< TensorMap33d >( stateVars, "Fp" );
    memcpy( Fp.data(), Spatial3D::I.data(), 9 * sizeof( double ) );
  }

  std::tuple< autodiff::dual, autodiff::dual, autodiff::dual > FiniteStrainGradientVonMises::fischerBurmeisterFunction(
    const autodiff::dual a,
    const autodiff::dual b,
    const double         epsilon ) const
  {
    const autodiff::dual sqrtTerm = sqrt( a * a + b * b + epsilon );
    const autodiff::dual f        = sqrtTerm - ( a + b );
    const autodiff::dual df_da    = 0.5 * a / sqrtTerm - 1.0;
    const autodiff::dual df_db    = 0.5 * b / sqrtTerm - 1.0;
    return { f, df_da, df_db };
  }

  autodiff::dual FiniteStrainGradientVonMises::fy( const autodiff::dual& kappa,
                                                   const autodiff::dual& laplaceKappa ) const
  {
    return fy0 + H * kappa - g * laplaceKappa;
  }

  void FiniteStrainGradientVonMises::computeStressAD( responseAD& res, const incrementAD& inc ) const
  {
    TensorMap33d                FpOld_map = stateLayout.getAs< TensorMap33d >( res.stateVars, "Fp" );
    Tensor33t< autodiff::dual > FpOld     = Marmot::makeDual( Tensor33d( FpOld_map ) );

    const autodiff::dual dLambda = inc.dLambda( 0 );

    // Accumulate kappa and laplaceKappa from state variables
    double& kappaOld        = stateLayout.getAs< double& >( res.stateVars, "kappa" );
    double& laplaceKappaOld = stateLayout.getAs< double& >( res.stateVars, "laplaceKappa" );

    autodiff::dual kappa        = kappaOld + inc.dLambda( 0 );
    autodiff::dual laplaceKappa = laplaceKappaOld + inc.laplaceDLambda( 0 );

    // 1. Trial state: F^{e,trial} = F * (F^{p,old})^{-1}
    Tensor33t< autodiff::dual > FeTrial = inc.deformation.F % Fastor::inverse( FpOld );
    Tensor33t< autodiff::dual > CeTrial = transpose( FeTrial ) % FeTrial;

    auto [psiTrial, dPsi_dCeTrial]     = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( CeTrial, K, G );
    Tensor33t< autodiff::dual > STrial = multiplyFastorTensorWithScalar( dPsi_dCeTrial, autodiff::dual( 2.0 ) );

    autodiff::dual              J        = determinant( inc.deformation.F );
    autodiff::dual              Jinv     = autodiff::dual( 1.0 ) / J;
    Tensor33t< autodiff::dual > tauTrial = FeTrial % STrial % transpose( FeTrial );
    Tensor33t< autodiff::dual > tTrial   = multiplyFastorTensorWithScalar( tauTrial, Jinv );

    // Deviatoric trial Cauchy stress and J2
    Tensor33t< autodiff::dual > tTrial_dev = Marmot::deviatoric( tTrial );
    autodiff::dual              J2_tTrial  = autodiff::dual( 0.5 ) * einsum_ij_ij_hardcoded( tTrial_dev, tTrial_dev );
    autodiff::dual              f_tr       = sqrt( autodiff::dual( 3.0 ) * J2_tTrial ) - fy( kappa, laplaceKappa );

    // Flow direction: df_p/dt = (3/2) * t_dev / sqrt(3 J2)
    // df_p/dM_{KL} = (df_p/dt_{ij}) * (1/J) * F^{e,-1}_{Ki} * F^{e}_{jL}
    Tensor33t< autodiff::dual > dfp_dt;
    if ( J2_tTrial > 1e-12 ) {
      autodiff::dual factor = autodiff::dual( 1.5 ) / sqrt( autodiff::dual( 3.0 ) * J2_tTrial );
      dfp_dt                = multiplyFastorTensorWithScalar( tTrial_dev, factor );
    }
    else {
      dfp_dt = Marmot::makeDual( Tensor33d( 0.0 ) );
    }

    Tensor33t< autodiff::dual > FeTrial_inv     = Fastor::inverse( FeTrial );
    Tensor33t< autodiff::dual > dfp_dM_unscaled = Tensor33t< autodiff::dual >( transpose( FeTrial_inv ) % dfp_dt %
                                                                               transpose( FeTrial ) );
    Tensor33t< autodiff::dual > dfp_dM          = multiplyFastorTensorWithScalar( dfp_dM_unscaled, Jinv );

    // 2. Exponential map: delta Fp = exp( dLambda * dfp_dM )
    Tensor33t< autodiff::dual > dGp     = multiplyFastorTensorWithScalar( dfp_dM, dLambda );
    Tensor33t< autodiff::dual > deltaFp = ContinuumMechanics::FiniteStrain::Plasticity::FlowIntegration::exponentialMap(
      dGp );

    Tensor33t< autodiff::dual > FpNew = deltaFp % FpOld;

    // 3. Updated elastic state: F^e = F * (F^{p,new})^{-1}
    Tensor33t< autodiff::dual > Fe  = inc.deformation.F % Fastor::inverse( FpNew );
    Tensor33t< autodiff::dual > Ce  = transpose( Fe ) % Fe;
    auto [psi, dPsi_dCe]            = EnergyDensityFunctions::FirstOrderDerived::PenceGouPotentialB( Ce, K, G );
    Tensor33t< autodiff::dual > S   = multiplyFastorTensorWithScalar( dPsi_dCe, autodiff::dual( 2.0 ) );
    Tensor33t< autodiff::dual > tau = Fe % S % transpose( Fe );

    // 4. Fischer-Burmeister complementarity function
    double scale             = 1e6;
    auto [fFB, df_da, df_db] = fischerBurmeisterFunction( -f_tr, dLambda * scale, 1e-12 );

    // 5. Update state variables (primal values only)
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        res.stateVars[i * 3 + j] = double( FpNew( i, j ) );
      }
    }
    kappaOld        = double( kappa );
    laplaceKappaOld = double( laplaceKappa );

    // 6. Write response
    res.tau                  = tau;
    res.f( 0 )               = fFB;
    res.elasticEnergyDensity = psi;
    res.dissipation          = autodiff::dual( 0.0 );
  }
} // namespace Marmot::Materials
