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

    /**
     * @brief Smoothing of the Fischer-Burmeister complementarity function, relative to @ref fy0.
     *
     * The function is @f$ \varphi(a,b) = \sqrt{a^2 + b^2 + \varepsilon} - (a+b) @f$ with
     * @f$ a = -f @f$ in stress units and @f$ b = \lambda\,\mathrm{scale} @f$, and
     * @f$ \varepsilon @f$ sets the radius of the corner at @f$ a = b = 0 @f$. Since @f$ a @f$ is a
     * stress, @f$ \varepsilon @f$ has to scale with the square of a stress; a fixed absolute value
     * cannot be right for every material. It is therefore stored relative to the initial yield
     * strength, @f$ \varepsilon = (\texttt{fbSmoothingRelative}\; f_{y0})^2 @f$.
     *
     * Why it must not be made arbitrarily small: @f$ \partial\varphi/\partial a @f$ multiplies the
     * whole yield condition row of the algorithmic tangent. At a converged plastic point
     * @f$ a = 0 @f$, so the distance to the corner is @f$ b @f$, which is proportional to the load
     * increment. Once @f$ b @f$ falls below the residual the solver is working at, a
     * residual-sized perturbation of @f$ a @f$ swings
     * @f$ \partial\varphi/\partial a @f$ from @f$ -1 @f$ to nearly zero, the yield condition
     * decouples from the displacements, and Newton oscillates instead of converging -- and cutting
     * the increment makes it worse, because that shrinks @f$ b @f$ further. Keeping
     * @f$ \sqrt{\varepsilon} @f$ above the working residual bounds that swing.
     */
    double fbSmoothingRelative = 1e-5;

    std::tuple< double, double, double > fy( double kappa, double laplaceKappa ) const;

    std::tuple< double, Vector6d, Matrix6d, double, double > yieldFunction( const Vector6d& stress,
                                                                            const double&   kappa,
                                                                            const double&   laplaceKappa ) const;

    std::tuple< double, double, double > fischerBurmeisterFunction( const double a,
                                                                    const double b,
                                                                    const double epsilon ) const;
  };

} // namespace Marmot::Materials
