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
#include "Marmot/MarmotFastorTensorBasics.h"

namespace Marmot {
  namespace ContinuumMechanics::FiniteStrain::Viscoelasticity {

    /**
     * @brief Struct to hold properties of a generalized maxwell model
     */
    struct MaxwellProperties {
      size_t                nMaxwell; ///< Number of maxwell elements
      std::vector< double > gamma;    ///< Weighting factors of maxwell elements
      double                sumGamma; ///< Sum of weighting factors
      std::vector< double > tau;      ///< Relaxation times of maxwell elements
    };

    MaxwellProperties createMaxwellProperties( int nMaxwell, const double* gammaTauPairVector );

    /**
     * @brief Evaluate generalized maxwell model contribution to stress and tangent
     * @param[in,out] stress             Stress tensor to be updated (input is the instantaneous hyperelastic stress)
     * @param[in,out] tangent            Tangent tensor to be updated (input is the instantaneous hyperelasticelastic
     * tangent)
     * @param[in]  dTangent_dDeformation Derivative of tangent w.r.t. deformation
     * @param[in]  initialCompliance     Initial compliance tensor
     * @param[in]  dStress               Increment of the initial stress tensor
     * @param[in]  dT                    Time increment
     * @param[in]  maxwellProperties     Properties of the generalized maxwell model
     * @param[in,out] stateVars          State variables array (size: nMaxwell * 9)
     *
     * @details This update function uses the approach in:
     * Liu et al. 2021, A continuum and computational framework for viscoelastodynamics: I. Finite deformation linear
     * models, Computer Methods Appl. Mech. Engrg. 385, 114059
     */
    void evaluateGeneralizedMaxwellModel( FastorStandardTensors::Tensor33d&           stress,
                                          FastorStandardTensors::Tensor3333d&         tangent,
                                          const FastorStandardTensors::Tensor333333d& dTangent_dDeformation,
                                          const FastorStandardTensors::Tensor3333d&   initialCompliance,
                                          const FastorStandardTensors::Tensor33d&     dStress,
                                          const double                                dT,
                                          const MaxwellProperties&                    maxwellProperties,
                                          double*                                     stateVars );
  } // namespace ContinuumMechanics::FiniteStrain::Viscoelasticity
} // namespace Marmot
