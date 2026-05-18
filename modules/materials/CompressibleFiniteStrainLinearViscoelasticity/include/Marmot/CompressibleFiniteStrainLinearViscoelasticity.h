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
#include <map>

namespace Marmot::Materials {

  /**
   * @class Marmot::Materials::CompressibleFiniteStrainLinearViscoelasticity
   * @brief Linear Viscoelastic Compressible Neo-Hookean hyperelastic material model
   * using the Pence–Gou potential, variant B.
   */
  class CompressibleFiniteStrainLinearViscoelasticity : public MarmotMaterialFiniteStrain {
  public:
    using MarmotMaterialFiniteStrain::MarmotMaterialFiniteStrain;

    /**
     * @brief Construct a CompressibleFiniteStrainLinearViscoelasticity material.
     * @param materialProperties Pointer to the material properties vector.
     * @param nMaterialProperties Length of @c materialProperties.
     * @param materialLabel Material label.
     */
    CompressibleFiniteStrainLinearViscoelasticity( const double* materialProperties,
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

    double getDensity( const double* stateVars ) const override;

  protected:
    /// @brief Enum for hyperelastic base
    enum HyperelasticBase { NeoHooke = 0, Yeoh = 1, MooneyRivlin = 2, PenceGouNeoHooke = 3 };

    /// map for hyperelastic base to number of elastic properties
    const std::map< HyperelasticBase, int > nElasticPropertiesMap = {
      { NeoHooke, 2 },
      { Yeoh, 4 },
      { MooneyRivlin, 3 },
      { PenceGouNeoHooke, 2 },
    };

    /// @brief Hyperelastic base model
    const HyperelasticBase hyperelasticBase;

    /// @brief shear only flag
    const double onlyShearCreep;

    /// @brief Elastic properties vector
    const Eigen::Map< const Eigen::VectorXd > elasticProperties;

    /// @brief Number of Maxwell elements
    const ContinuumMechanics::FiniteStrain::Viscoelasticity::MaxwellProperties maxwellProperties;

    /// @brief Initial compliance tensor
    FastorStandardTensors::Tensor3333d initialCompliance;

    void initializeStateLayout()
    {
      stateLayout.add( "S0_old", 9 );
      stateLayout.add( "creepStateVars", maxwellProperties.nMaxwell * 9 ); // plastic deformation gradient
      stateLayout.finalize();
    }

    std::tuple< double,
                FastorStandardTensors::Tensor33d,
                FastorStandardTensors::Tensor3333d,
                FastorStandardTensors::Tensor333333d >
    computeEnergyDensityAndDerivatives( const FastorStandardTensors::Tensor33d& C ) const;
  };

} // namespace Marmot::Materials
