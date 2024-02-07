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
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/WMPPlasticity.h"
#include <iostream>
#include <string>
#include <vector>

namespace Marmot::Materials {
  /**
   * \brief Implementation of the multisurface plasticity model by Schmidt & Kaliske (2006).
   *
   */
  class WoodMultisurfacePlasticity : public MarmotMaterialHypoElastic {
  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

    WoodMultisurfacePlasticity( const double* materialProperties, int nMaterialProperties, int materialNumber );

  protected:
    Matrix6d transformationMatrixStrain, transformationMatrixStrainInv, transformationMatrixStress,
      transformationMatrixStressInv, localElasticStiffnessTensor, localElasticComplianceTensor,
      globalElasticStiffnessTensor;
    Matrix3d localCoordinateSystem;

    // elastic constants
    const double &ER, ET, EL;
    const double &nuTR, nuLR, nuLT;
    const double &GRT, GRL, GTL;

    // material coordinate system
    const Vector3d nR, nT;

    // plasticity
    WMPPlasticity plasticity;

    void computeStress( double* stress,
                        double* dStressDDStrain,

                        const double* dStrain,
                        const double* timeOld,
                        const double  dT,
                        double&       pNewDT ) override;

    class WMPStateVarManager : public MarmotStateVarVectorManager {

    public:
      inline const static auto layout = makeLayout( {
        { .name = "alpha", .length = 7 },
      } );

      Eigen::Map< Eigen::Vector< double, 7 > > alpha;

      WMPStateVarManager( double* theStateVarVector )
        : MarmotStateVarVectorManager( theStateVarVector, layout ), alpha( &find( "alpha" ) ){};
    };
    std::unique_ptr< WMPStateVarManager > managedStateVars;

    StateView getStateView( const std::string& result ) override;

    int getNumberOfRequiredStateVars() override { return WMPStateVarManager::layout.nRequiredStateVars; }

    void assignStateVars( double* stateVars, int nStateVars ) override;
  };
} // namespace Marmot::Materials
