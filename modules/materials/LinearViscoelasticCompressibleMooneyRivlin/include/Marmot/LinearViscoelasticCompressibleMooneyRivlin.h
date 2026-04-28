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
 * Alexander Dummer alexander.dummer@uibk.ac.at
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
#include "Marmot/MarmotFiniteStrainViscoelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"

namespace Marmot::Materials {

  /**
   * @class Marmot::Materials::LinearViscoelasticCompressibleMooneyRivlin
   * @brief Linear viscoelastic compressible Mooney–Rivlin hyperelastic material model
   * using a Mooney–Rivlin type strain-energy potential.
   */
  class LinearViscoelasticCompressibleMooneyRivlin : public MarmotMaterialFiniteStrain {
  public:
    using MarmotMaterialFiniteStrain::MarmotMaterialFiniteStrain;

    /**
     * @brief Construct a LinearViscoelasticCompressibleMooneyRivlin material.
     * @param materialProperties Pointer to the material properties vector.
     * @param nMaterialProperties Length of @c materialProperties.
     * @param materialLabel Material label.
     */
    LinearViscoelasticCompressibleMooneyRivlin( const double* materialProperties,
                                                int           nMaterialProperties,
                                                int           materialLabel );

    /**
     * @brief Compute the Kirchhoff stress and the algorithmic tangent for the current step.
     *
     * @param[in,out] response
     *   - @c tau - Kirchhoff stress tensor @f$\boldsymbol{\tau}@f$.
     *   - @c elasticEnergyDensity - elastic energy density  @f$\psi@f$.
     *   - @c rho - density (unused here).
     * @param[in,out] tangents
     *   - @c dTau_dF - algorithmic tangent @f$\partial\boldsymbol{\tau}/\partial\boldsymbol{F}@f$.
     * @param[in]  deformation
     *   - @c F - deformation gradient @f$\boldsymbol{F}@f$.
     * @param[in]  timeIncrement
     *   - @c t - old time.
     *   - @c dT - time increment.
     *
     * Template parameter @c <3> indicates 3D.
     */
    void computeStress( ConstitutiveResponse< 3 >&,
                        AlgorithmicModuli< 3 >&,
                        const Deformation< 3 >&,
                        const TimeIncrement& ) const override;

  protected:
    /// @brief Mooney-Rivlin material parameter C1
    const double& C1;
    /// @brief Mooney-Rivlin material parameter C2
    const double& C2;
    /// @brief Bulk modulus
    const double& K;
    /// @brief Number of Maxwell elements
    const ContinuumMechanics::FiniteStrain::Viscoelasticity::MaxwellProperties maxwellProperties;
    /// @brief Initial compliance tensor
    FastorStandardTensors::Tensor3333d initialCompliance;
    void                               initializeStateLayout() override
    {
      stateLayout.add( "S0_old", 9 );
      stateLayout.add( "creepStateVars", maxwellProperties.nMaxwell * 9 ); // plastic deformation gradient
      stateLayout.finalize();
    }
  };

} // namespace Marmot::Materials
