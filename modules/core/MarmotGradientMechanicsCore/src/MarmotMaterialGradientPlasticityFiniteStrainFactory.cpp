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
#include "Marmot/MarmotMaterialGradientPlasticityFiniteStrainFactory.h"
#include <iostream>
#include <stdexcept>

namespace MarmotLibrary {

  template < int nYieldSurfaces >
  MarmotMaterialGradientPlasticityFiniteStrain< nYieldSurfaces >* MarmotMaterialGradientPlasticityFiniteStrainFactory<
    nYieldSurfaces >::createMaterial( const std::string& materialName,
                                      const double*      materialProperties,
                                      int                nMaterialProperties,
                                      int                materialNumber )
  {
    auto& map = materialFactoryFunctionByName();

    auto it = map.find( materialName );
    if ( it != map.end() ) {
      return it->second( materialProperties, nMaterialProperties, materialNumber );
    }
    else {
      throw std::invalid_argument( "Material " + materialName + " not registered or found!" );
    }
  }

  template < int nYieldSurfaces >
  typename MarmotMaterialGradientPlasticityFiniteStrainFactory< nYieldSurfaces >::MaterialFactoryMap& MarmotMaterialGradientPlasticityFiniteStrainFactory<
    nYieldSurfaces >::materialFactoryFunctionByName()
  {
    static MaterialFactoryMap map;
    return map;
  }

  // Explicit instantiation for common cases
  template class MarmotMaterialGradientPlasticityFiniteStrainFactory< 1 >;
  template class MarmotMaterialGradientPlasticityFiniteStrainFactory< 2 >;
  template class MarmotMaterialGradientPlasticityFiniteStrainFactory< 3 >;

} // namespace MarmotLibrary
