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
#include "Marmot/MarmotMaterialFiniteStrainAD.h"

namespace Marmot::Materials {

  /**
   * @class Marmot::Materials::FiniteStrainOrthotropicBiotViscoelasticity
   * @brief Linear Viscoelastic Biot hyperelastic material model
   * considering orthotropic elasticity and linear viscoelasticity in the Biot stress, i.e. the stress is computed as
   */
  class FiniteStrainOrthotropicBiotViscoelasticity : public MarmotMaterialFiniteStrainAD {
  public:
    using MarmotMaterialFiniteStrainAD::MarmotMaterialFiniteStrainAD;

    /**
     * @brief Construct a FiniteStrainOrthotropicBiotViscoelasticity material.
     * @param materialProperties Pointer to the material properties vector.
     * @param nMaterialProperties Length of @c materialProperties.
     * @param materialLabel Material label.
     */
    FiniteStrainOrthotropicBiotViscoelasticity( const double* materialProperties,
                                                int           nMaterialProperties,
                                                int           materialLabel );

    void computeStressAD( ConstitutiveResponseAD< 3 >&,
                          const DeformationAD< 3 >&,
                          const TimeIncrement& ) const override;

    double getDensity( const double* stateVars ) const override;

  protected:
    const double E1;
    const double E2;
    const double E3;

    const double nu12;
    const double nu13;
    const double nu23;

    const double G12;
    const double G13;
    const double G23;

    const ContinuumMechanics::FiniteStrain::Viscoelasticity::MaxwellProperties maxwellProperties;

    FastorStandardTensors::Tensor3333d dBiotStress_dU = FastorStandardTensors::Tensor3333d( 0 );

    void initializeStateLayout()
    {
      stateLayout.add( "S0_old", 9 );
      stateLayout.add( "creepStateVars", maxwellProperties.nMaxwell * 9 ); // plastic deformation gradient
      stateLayout.finalize();
    }
  };

} // namespace Marmot::Materials
