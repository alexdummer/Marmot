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
 * Modified for C0-Continuous Penalty-Enhanced Gradient Plasticity
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

#include "Marmot/MarmotMaterialC0PenaltyGradientPlasticityFactory.h"
#include "Marmot/MarmotJournal.h"

using namespace MarmotLibrary;

MarmotMaterialC0PenaltyGradientPlasticityFactory::MaterialFactoryMap& MarmotMaterialC0PenaltyGradientPlasticityFactory::
  materialFactoryFunctionByName()
{
  static MaterialFactoryMap map;
  return map;
}

MarmotMaterialC0PenaltyGradientPlasticity* MarmotMaterialC0PenaltyGradientPlasticityFactory::createMaterial(
  const std::string& materialName,
  const double*      materialProperties,
  int                nMaterialProperties,
  int                materialNumber )
{
  auto& map = materialFactoryFunctionByName();
  auto  it  = map.find( materialName );
  if ( it == map.end() ) {
    std::string reg = "Registered materials are: ";
    for ( const auto& pair : map ) {
      reg += pair.first + ", ";
    }
    throw std::invalid_argument( MakeString()
                                 << __PRETTY_FUNCTION__ << " Material " + materialName + " not registered!" + reg );
  }

  return it->second( materialProperties, nMaterialProperties, materialNumber );
}
