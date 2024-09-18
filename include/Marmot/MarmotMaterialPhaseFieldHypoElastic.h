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
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * Magdalena Schreter magdalena.schreter@uibk.ac.at
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
#include "Marmot/MarmotMaterialPhaseFieldMechanical.h"

class MarmotMaterialPhaseFieldHypoElastic : public MarmotMaterialPhaseFieldMechanical {

public:
  using MarmotMaterialPhaseFieldMechanical::MarmotMaterialPhaseFieldMechanical;

  virtual void computeStress( double*       stress,
                              double*       KLocal,
                              double&       nonLocalRadius,
                              double*       dStress_dDeformationGradient,
                              double*       dKLocal_dDeformationGradient,
                              double*       dStress_dK,
                              double*       dKlocal_dK,
                              const double* FOld,
                              const double* FNew,
                              const double* KOld,
                              const double* dK,
                              const double time,
                              const double  dT) override;

  virtual void computeStress( double*       stress,
                              double*       K_local,
                              double&       nonLocalRadius,
                              double*       dStressDDStrain,
                              double*       dK_localDDStrain,
                              double*       dStressDK,
                              double*       dKlocal_dK,
                              const double* dStrain,
                              const double*       KOld,
                              const double*       dK,
                              const double time,
                              const double  dT) = 0;

  /* using MarmotMaterialPhaseFieldMechanical::computePlaneStress; */
  /* virtual void computePlaneStress( double*       stress2D, */
  /*                                  double*       KLocal2D, */
  /*                                  double&       nonLocalRadius, */
  /*                                  double*       dStress_DStrain2D, */
  /*                                  double*       dKLocal_dStrain2D, */
  /*                                  double*       dStress_dK2D, */
  /*                                  const double* dStrain2D, */
  /*                                  const double*       KOld, */
  /*                                  const double*       dK, */
  /*                                  const double time, */
  /*                                  const double  dT); */
};
