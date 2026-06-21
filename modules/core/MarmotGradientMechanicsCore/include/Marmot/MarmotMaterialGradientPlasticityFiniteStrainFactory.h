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
#include "Marmot/MarmotMaterialGradientPlasticityFiniteStrain.h"
#include <cassert>
#include <functional>
#include <string>
#include <unordered_map>

namespace MarmotLibrary {

  /**
   * @brief Factory class for creating MarmotMaterialGradientPlasticityFiniteStrain instances.
   */
  template < int nYieldSurfaces >
  class MarmotMaterialGradientPlasticityFiniteStrainFactory {
  public:
    using materialFactoryFunction = std::function< MarmotMaterialGradientPlasticityFiniteStrain<
      nYieldSurfaces >*( const double* materialProperties, int nMaterialProperties, int materialNumber ) >;

    MarmotMaterialGradientPlasticityFiniteStrainFactory() = delete;

    static MarmotMaterialGradientPlasticityFiniteStrain< nYieldSurfaces >* createMaterial(
      const std::string& materialName,
      const double*      materialProperties,
      int                nMaterialProperties,
      int                materialNumber );

    template < class T >
    static bool registerMaterial( const std::string& materialName )
    {
      auto& map = materialFactoryFunctionByName();

      assert( map.find( materialName ) == map.end() && "Material already registered!" );

      map[materialName] = []( const double* materialProperties,
                              int           nMaterialProperties,
                              int materialNumber ) -> MarmotMaterialGradientPlasticityFiniteStrain< nYieldSurfaces >* {
        return new T( materialProperties, nMaterialProperties, materialNumber );
      };
      return true;
    }

    using MaterialFactoryMap = std::unordered_map< std::string, materialFactoryFunction >;
    static MaterialFactoryMap& materialFactoryFunctionByName();
  };

} // namespace MarmotLibrary

#define REGISTER_MARMOT_MATERIAL_GRADIENT_FINITE_STRAIN( MATERIAL_CLASS_NAME, YIELD_SURFACES )                         \
  bool is##MATERIAL_CLASS_NAME##Registered = MarmotLibrary::MarmotMaterialGradientPlasticityFiniteStrainFactory<       \
    YIELD_SURFACES >::registerMaterial< MATERIAL_CLASS_NAME >( #MATERIAL_CLASS_NAME );
