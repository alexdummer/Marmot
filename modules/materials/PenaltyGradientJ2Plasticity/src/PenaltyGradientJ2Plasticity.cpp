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
 * Modified for C0-Continuous Penalty-Enhanced Gradient Plasticity
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

#include "Marmot/PenaltyGradientJ2Plasticity.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;

  namespace {
    constexpr double innerNewtonTol        = 1e-12;
    constexpr int    nMaxInnerNewtonCycles = 15;
  } // namespace

  PenaltyGradientJ2Plasticity::PenaltyGradientJ2Plasticity( const double* materialProperties,
                                                            int           nMaterialProperties,
                                                            int           materialNumber )
    : MarmotMaterialC0PenaltyGradientPlasticity( materialProperties, nMaterialProperties, materialNumber ),
      Cel( ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( materialProperties[0], materialProperties[1] ) )
  {
    stateLayout.add( "kappa", 1 );
    stateLayout.finalize();
  }

  double PenaltyGradientJ2Plasticity::getDensity( const double* stateVars ) const
  {
    if ( nMaterialProperties < 7 ) {
      throw std::runtime_error(
        std::string( MakeString() << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 7." ) );
    }
    return materialProperties[6];
  }

  void PenaltyGradientJ2Plasticity::computeStress( response& res, tangents& tan, const increment& inc ) const
  {
    // Material properties
    const double& E                = materialProperties[0];
    const double& nu               = materialProperties[1];
    const double& yieldStress      = materialProperties[2];
    const double& HLin             = materialProperties[3];
    const double& deltaYieldStress = materialProperties[4];
    const double& delta            = materialProperties[5];

    // Shear modulus
    const double G = E / ( 2. * ( 1. + nu ) );

    // Access state variables: local cumulative plastic strain
    double& kappa = stateLayout.getAs< double& >( res.stateVars, "kappa" );

    // Isotropic hardening law using nonlocal cumulative plastic strain
    auto fy = [&]( double kappaNL_ ) {
      return yieldStress + HLin * kappaNL_ + deltaYieldStress * ( 1. - std::exp( -delta * kappaNL_ ) );
    };

    // Derivative of fy w.r.t. nonlocal kappa
    auto dfy_dKappaNL = [&]( double kappaNL_ ) {
      return HLin + deltaYieldStress * delta * std::exp( -delta * kappaNL_ );
    };

    // Map stress
    mVector6d S( res.stress.data() );

    // Compute elastic predictor
    const Vector6d trialStress = S + Cel * inc.dStrain;

    using namespace ContinuumMechanics::VoigtNotation;
    const double rhoTrial = std::sqrt( 2. * Invariants::J2( trialStress ) );

    // Yield function: f = rho - sqrt(2/3) * fy(kappaNL)
    const double fTrial = rhoTrial - Constants::sqrt2_3 * fy( inc.kappaNL );

    if ( fTrial >= 0.0 ) {
      // Plastic step - return mapping

      auto g = [&]( double deltaKappa ) {
        return rhoTrial - Constants::sqrt6 * G * deltaKappa - Constants::sqrt2_3 * fy( inc.kappaNL );
      };

      // Return mapping direction
      Vector6d n = ContinuumMechanics::VoigtNotation::IDev * trialStress / rhoTrial;

      // Since the yield function depends on kappaNL (not the local kappa),
      // the return mapping for dKappa is direct:
      //   g(dKappa) = rhoTrial - sqrt(6)*G*dKappa - sqrt(2/3)*fy(kappaNL) = 0
      //   dKappa = (rhoTrial - sqrt(2/3)*fy(kappaNL)) / (sqrt(6)*G)
      double dKappa = ( rhoTrial - Constants::sqrt2_3 * fy( inc.kappaNL ) ) / ( Constants::sqrt6 * G );

      if ( dKappa < 0. )
        dKappa = 0.;

      double dLambda = Constants::sqrt3_2 * dKappa;

      // Update stress
      S = trialStress - 2. * G * dLambda * n;

      // Update local cumulative plastic strain
      kappa += dKappa;

      // Set local kappa in response
      res.kappaLocal = kappa;

      // --- Algorithmic tangents ---

      // Consistent tangent dStress/dStrain
      Matrix6d IDevHalfShear = ContinuumMechanics::VoigtNotation::IDev;
      IDevHalfShear.block< 6, 3 >( 0, 3 ) *= 0.5;

      tan.dStressDDStrain = Cel - 4. * G * G * dLambda / rhoTrial * IDevHalfShear;

      // dStress/dKappaNL: stress depends on kappaNL through dLambda
      // dLambda/dKappaNL = sqrt(3/2) * dKappa/dKappaNL
      // dKappa/dKappaNL = -sqrt(2/3) * dfy/dKappaNL / (sqrt(6)*G)
      //                 = -dfy/dKappaNL / (3*G)
      const double dKappa_dKappaNL  = -Constants::sqrt2_3 * dfy_dKappaNL( inc.kappaNL ) / ( Constants::sqrt6 * G );
      const double dLambda_dKappaNL = Constants::sqrt3_2 * dKappa_dKappaNL;

      tan.dStressDDKappaNL = -2. * G * dLambda_dKappaNL * n;

      // dKappaLocal/dStrain: dKappa/dStrain
      // dKappa = (rhoTrial - sqrt(2/3)*fy(kappaNL)) / (sqrt(6)*G)
      // dKappa/dStrain = (1/(sqrt(6)*G)) * dRhoTrial/dStrain
      // dRhoTrial/dStrain = (1/rhoTrial) * s^T * Cel (where s is deviatoric trial stress)
      Vector6d dRho_dStrain = ( 1.0 / rhoTrial ) *
                              ( ContinuumMechanics::VoigtNotation::IDev * trialStress ).transpose() * Cel;
      // Correct for Voigt: need to account for factor of 2 in shear terms of J2
      // J2 = 0.5 * s^T * P * s where P accounts for Voigt factors
      // Actually: rhoTrial = sqrt(2*J2), dRho/dStrain = (1/rhoTrial) * dJ2/d(eps_trial) * d(eps_trial)/d(dStrain)
      // But eps_trial is the elastic strain, so d(trial_stress)/d(dStrain) = Cel

      // Simplified: using the chain rule with Voigt conventions
      Vector6d sVoigt       = ContinuumMechanics::VoigtNotation::IDev * trialStress;
      Vector6d sVoigtScaled = sVoigt;
      sVoigtScaled.tail( 3 ) *= 2.0; // Account for Voigt convention in J2 derivative

      tan.dKappaLocalDDStrain = ( 1. / ( Constants::sqrt6 * G * rhoTrial ) ) * sVoigtScaled.transpose() * Cel;

      // dKappaLocal/dKappaNL
      tan.dKappaLocalDDKappaNL = dKappa_dKappaNL;
    }
    else {
      // Elastic step
      S = trialStress;

      // Local kappa unchanged
      res.kappaLocal = kappa;

      // Tangents
      tan.dStressDDStrain = Cel;
      tan.dStressDDKappaNL.setZero();
      tan.dKappaLocalDDStrain.setZero();
      tan.dKappaLocalDDKappaNL = 0.0;
    }
  }

} // namespace Marmot::Materials
