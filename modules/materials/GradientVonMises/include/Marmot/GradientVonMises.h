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

#pragma once
#include "Marmot/MarmotMaterialGradientPlasticityHypoElastic.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

namespace Marmot::Materials {

  class GradientVonMises : public MarmotMaterialGradientPlasticityHypoElastic< 1 > {

  public:
    GradientVonMises( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void computeStress( response& res, tangents& tan, const increment& inc ) const override;

    void initializeStateLayout()
    {
      stateLayout.add( "kappa", 1 );
      stateLayout.add( "laplaceKappa", 1 );
      stateLayout.finalize();
    }

    double getDensity( const double* stateVars ) const override;

    std::vector< double > getNonlocalViscosity( const double* stateVars ) const override;

  private:
    /// @brief Elastic stiffness tensor
    const Marmot::Matrix6d C;
    const double&          fy0; // initial yield strength
    const double&          H;   // hardening modulus
    const double&          g;   // gradient influence parameter

    /// compute yield stress and its derivatives from kappa and laplaceKappa
    std::tuple< double, double, double > fy( double kappa, double laplaceKappa ) const
    {
      const double sigmaY = fy0 + H * kappa - g * laplaceKappa; // yield stress as a function of kappa
      return { sigmaY,
               H,
               -g }; // return yield stress and its derivatives (d(sigmaY)/d(kappa) = H, d(sigmaY)/d(laplaceKappa) = -g)
    }

    // compute the von Mises yield function value and its derivatives with respect to stress, kappa, and laplaceKappa
    std::tuple< double, Vector6d, Matrix6d, double, double > yieldFunction( const Vector6d& stress,
                                                                            double&         kappa,
                                                                            double&         laplaceKappa ) const
    {
      using namespace Marmot::ContinuumMechanics::VoigtNotation;
      const auto [sigmaY, dSigmaY_dKappa, dSigmaY_dLaplaceKappa] = fy( kappa, laplaceKappa );
      const double   J2                                          = Invariants::J2( stress );
      const Vector6d dJ2_dStress                                 = Derivatives::dJ2_dStress( stress );
      const Matrix6d d2J2_dStress2                               = Derivatives::d2J2_dStress2( stress );

      const Vector6d dF_ddStress  = dJ2_dStress * ( 3.0 / ( 2.0 * std::sqrt( 3.0 * J2 ) + 1e-14 ) );
      const Matrix6d d2F_dStress2 = d2J2_dStress2 * ( 3.0 / ( 2.0 * std::sqrt( 3.0 * J2 ) + 1e-14 ) ) -
                                    ( dJ2_dStress * dJ2_dStress.transpose() ) * ( 9.0 / ( 8.0 * std::pow( J2, 1.5 ) ) );

      const double f = std::sqrt( 3.0 * J2 ) - sigmaY; // yield

      return { f,
               dF_ddStress,
               d2F_dStress2,
               -dSigmaY_dKappa,
               -dSigmaY_dLaplaceKappa }; // return yield function value and its derivatives
    }
  };

} // namespace Marmot::Materials
