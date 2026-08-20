/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
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
#include "Marmot/MarmotMaterialFiniteStrainAD.h"

namespace Marmot::Materials {

  /**
   * @class Marmot::Materials::IncompressibleNeoHooke
   * @brief Classical incompressible (purely isochoric) Neo-Hookean hyperelastic material model.
   *
   * @par Material parameters
   * - @b mu - shear modulus (materialProperties[0])
   * - @b density (optional) - materialProperties[1]
   *
   * @par State variables
   * - No state variables required.
   *
   * @details Implements the standard isochoric Neo-Hookean potential
   * \f[
   *   \Psi_\mathrm{iso}(\bar\lambda) = \frac{\mu}{2} \left( \bar\lambda_1^2 + \bar\lambda_2^2 +
   *   \bar\lambda_3^2 - 3 \right) = \frac{\mu}{2} \left( \bar I_1 - 3 \right),
   * \f]
   * with isochoric principal stretches \f$\bar\lambda_i = J^{-1/3}\lambda_i\f$ and isochoric first
   * invariant \f$\bar I_1 = \bar\lambda_1^2+\bar\lambda_2^2+\bar\lambda_3^2\f$. This is exactly the
   * one-term Ogden model with exponent 2, and is implemented as such (see
   * Marmot::Materials::Ogden), reusing the same validated spectral machinery instead of
   * duplicating it.
   *
   * There is no volumetric energy term \f$U(J)\f$: this material resists no volumetric deformation
   * on its own and must therefore not be used with a plain displacement element. It is intended
   * exclusively for use with a mixed displacement-pressure element that enforces incompressibility
   * as an independent constraint, e.g. Marmot::Elements::DisplacementPressureFiniteStrainElement.
   *
   * @ingroup materials_hyperelastic
   */
  class IncompressibleNeoHooke : public MarmotMaterialFiniteStrainAD {

  public:
    /**
     * @brief Construct an IncompressibleNeoHooke material.
     * @param materialProperties Expects @c mu at index 0.
     * @param nMaterialProperties Length of @c materialProperties.
     * @param materialLabel Material label.
     */
    IncompressibleNeoHooke( const double* materialProperties, int nMaterialProperties, int materialLabel );

    /**
     * @brief Compute the isochoric Kirchhoff stress with dual numbers.
     * @param[in,out] response
     *   - @c tau - Kirchhoff stress tensor @f$\boldsymbol{\tau}@f$ (trace-free/deviatoric).
     *   - @c elasticEnergyDensity - isochoric elastic energy density @f$\Psi_\mathrm{iso}@f$.
     * @param[in]  deformation
     *   - @c F - deformation gradient @f$\boldsymbol{F}@f$.
     * @param[in]  timeIncrement Unused (hyperelastic, no history).
     */
    void computeStressAD( ConstitutiveResponseAD< 3 >& response,
                          const DeformationAD< 3 >&    deformation,
                          const TimeIncrement&         timeIncrement ) const override;

    double getDensity( const double* stateVars ) const override
    {
      if ( this->nMaterialProperties < 2 )
        throw std::runtime_error(
          std::string( MakeString() << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 2." ) );
      return this->materialProperties[1];
    }
  };

} // namespace Marmot::Materials
