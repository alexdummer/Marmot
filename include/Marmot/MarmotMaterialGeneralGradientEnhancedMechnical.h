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
#include "Marmot/MarmotMaterial.h"

class MarmotMaterialGeneralGradientEnhancedMechanical : public MarmotMaterial {

public:
  using MarmotMaterial::MarmotMaterial;

  virtual void computeStress( double*       stress,
                              double*       K_local,
                              double&       nonLocalRadius,
                              double*       dStressDDFNew,
                              double*       dK_localDDFNew,
                              double*       dStressDK,
                              double*       dKlocal_dK,
                              const double* FOld,
                              const double* FNew,
                              const double* KOld,
                              const double* dK,
                              const double  time,
                              const double  dT ) = 0;

  virtual void computePlaneStress( double*       stress,
                                   double*       K_local,
                                   double&       nonLocalRadius,
                                   double*       dStressDDFNew,
                                   double*       dK_localDDFNew,
                                   double*       dStressDK,
                                   double*       dKlocal_dK,
                                   const double* FOld,
                                   const double* FNew,
                                   double*       KOld,
                                   double*       dK,
                                   const double  time,
                                   const double  dT ) {};

  virtual void computeUniaxialStress( double*       stress,
                                      double*       K_local,
                                      double&       nonLocalRadius,
                                      double*       dStressDDFNew,
                                      double*       dK_localDDFNew,
                                      double*       dStressDK,
                                      const double* FOld,
                                      const double* FNew,
                                      double*       KOld,
                                      double*       dK,
                                      const double  time,
                                      const double  dT ) {};
};
