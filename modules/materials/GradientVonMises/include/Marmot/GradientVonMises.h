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
    void computeStressStandard( response& res, tangents& tan, const increment& inc ) const;
    void computeStressFischerBurmeister( response& res, tangents& tan, const increment& inc ) const;

    void initializeStateLayout()
    {
      stateLayout.add( "kappa", 1 );
      stateLayout.add( "laplaceKappa", 1 );
      stateLayout.finalize();
    }

    double getDensity( const double* stateVars ) const override;

    std::vector< double > getNonlocalViscosity( const double* stateVars ) const override;

  private:
    const double& E; // Young's modulus
    /// @brief Elastic stiffness tensor
    const Marmot::Matrix6d C;
    const double&          fy0; // initial yield strength
    const double&          H;   // hardening modulus
    const double&          g;   // gradient influence parameter

    enum class Implementation { standard, fischer_burmeister } implementation = Implementation::fischer_burmeister;

    std::tuple< double, double, double > fy( double kappa, double laplaceKappa ) const;

    std::tuple< double, Vector6d, Matrix6d, double, double > yieldFunction( const Vector6d& stress,
                                                                            const double&   kappa,
                                                                            const double&   laplaceKappa ) const;

    std::tuple< double, double, double > fischerBurmeisterFunction( const double a,
                                                                    const double b,
                                                                    const double epsilon ) const;
  };

} // namespace Marmot::Materials
