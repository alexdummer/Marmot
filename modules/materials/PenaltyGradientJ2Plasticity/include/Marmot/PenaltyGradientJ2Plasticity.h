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

#pragma once
#include "Marmot/MarmotMaterialC0PenaltyGradientPlasticity.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  /**
   * @brief Gradient-enhanced J2 (von Mises) plasticity for the C0-continuous
   *        penalty-enhanced gradient plasticity framework.
   *
   * @details Implements a J2 plasticity model with isotropic hardening (linear +
   * exponential saturation). The yield function depends on the nonlocal cumulative
   * plastic strain \f$ \bar{\kappa} \f$ rather than the local one, providing
   * mesh-independent results in softening regimes.
   *
   * The yield function reads:
   * \f[
   *   f = \| \mathbf{s} \| - \sqrt{\frac{2}{3}}\, \sigma_y(\bar{\kappa})
   * \f]
   * where \f$ \mathbf{s} \f$ is the deviatoric stress and
   * \f[
   *   \sigma_y(\bar{\kappa}) = \sigma_{y0} + H \bar{\kappa}
   *     + \Delta\sigma_y \left(1 - e^{-\delta \bar{\kappa}}\right)
   * \f]
   *
   * **Material properties** (in order):
   *  - \f$ E \f$             : Young's modulus
   *  - \f$ \nu \f$           : Poisson's ratio
   *  - \f$ \sigma_{y0} \f$   : Initial yield stress
   *  - \f$ H \f$             : Linear hardening modulus
   *  - \f$ \Delta\sigma_y \f$: Saturation stress increment
   *  - \f$ \delta \f$        : Saturation exponent
   *  - \f$ \rho \f$          : Mass density (optional, mandatory for dynamic simulations)
   */
  class PenaltyGradientJ2Plasticity : public MarmotMaterialC0PenaltyGradientPlasticity {

  public:
    PenaltyGradientJ2Plasticity( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void computeStress( response& res, tangents& tan, const increment& inc ) const override;

    double getDensity( const double* stateVars ) const override;

  private:
    /// @brief Elastic stiffness tensor
    const Marmot::Matrix6d Cel;
  };

} // namespace Marmot::Materials
