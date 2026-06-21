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

#include "Marmot/MarmotMaterialGradientPlasticityFiniteStrainAD.h"
#include "Marmot/MarmotPhaseFieldEnergyDegradation.h"

namespace Marmot::Materials {

  class FiniteStrainGradientVonMises : public MarmotMaterialGradientPlasticityFiniteStrainAD< 1 > {
  public:
    const double K;   ///< Bulk modulus
    const double G;   ///< Shear modulus
    const double fy0; ///< Initial yield stress
    const double H;   ///< Hardening modulus
    const double g;   ///< Gradient hardening modulus

    FiniteStrainGradientVonMises( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void computeStressAD( responseAD& res, const incrementAD& inc ) const override;

    double getDensity( const double* stateVars ) const override;

    std::vector< double > getNonlocalViscosity( const double* stateVars ) const override;

    void initializeYourself( double* stateVars, int nStateVars ) override;

  private:
    autodiff::dual fischerBurmeisterFunction( const autodiff::dual a,
                                              const autodiff::dual b,
                                              const double         epsilon ) const;

    autodiff::dual fy( const autodiff::dual& kappa, const autodiff::dual& laplaceKappa ) const;
  };

} // namespace Marmot::Materials
