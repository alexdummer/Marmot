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
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTypedefs.h"
#include <iostream>
#include <string>
#include <vector>

namespace Marmot::Materials {
  /**
   * \brief Implementation of the Model proposed by Hassani et al. (2015) restricted to constant humidity.
   *
   */
  class WoodCreepPlasticity : public MarmotMaterialHypoElastic {
  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

    WoodCreepPlasticity( const double* materialProperties, int nMaterialProperties, int materialNumber );

  protected:
    Matrix6d transformationMatrixStrain, transformationMatrixStress, localElasticStiffnessTensor;
    Matrix3d localCoordinateSystem;

    // elastic constants
    const double &ER, ET, EL;
    const double &nuTR, nuLR, nuLT;
    const double &GRT, GRL, GTL;

    // material coordinate system
    const Vector3d nR, nT;

    void computeStress( double* stress,
                        double* dStressDDStrain,

                        const double* dStrain,
                        const double* timeOld,
                        const double  dT,
                        double&       pNewDT );

    StateView getStateView( const std::string& result ) { return { nullptr, 0 }; };

    int getNumberOfRequiredStateVars() { return 0; }
  };
} // namespace Marmot::Materials
