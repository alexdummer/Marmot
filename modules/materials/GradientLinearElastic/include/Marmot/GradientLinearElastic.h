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

  class GradientLinearElastic : public MarmotMaterialGradientPlasticityHypoElastic< 1 > {

  public:
    GradientLinearElastic( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void computeStress( response& res, tangents& tan, const increment& inc ) const override;

    void initializeStateLayout() { stateLayout.finalize(); }

    double getDensity( const double* stateVars ) const override;

    std::vector< double > getNonlocalViscosity( const double* stateVars ) const override;

  private:
    /// @brief Young's modulus
    const double& E;
    /// @brief Poisson's ratio
    const double& nu;
    /// @brief First Lame parameter
    const double lambda;
    /// @brief Second Lame parameter
    const double mu;
  };

} // namespace Marmot::Materials
