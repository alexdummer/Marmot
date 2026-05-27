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
#include <cassert>
#include <functional>
#include <string>
#include <unordered_map>

namespace MarmotLibrary {

  class MarmotMaterialGradientPlasticityHypoElasticFactory {
  public:
    using materialFactoryFunction =
      std::function< Marmot::MarmotMaterialGradientPlasticityHypoElastic*( const double* materialProperties,
                                                                           int           nMaterialProperties,
                                                                           int           materialNumber ) >;

    MarmotMaterialGradientPlasticityHypoElasticFactory() = delete;

    static Marmot::MarmotMaterialGradientPlasticityHypoElastic* createMaterial( const std::string& materialName,
                                                                                 const double*      materialProperties,
                                                                                 int                nMaterialProperties,
                                                                                 int                materialNumber );

    template < class T >
    static bool registerMaterial( const std::string& materialName )
    {
      auto& map = materialFactoryFunctionByName();

      assert( map.find( materialName ) == map.end() && "Material already registered!" );

      map[materialName] = []( const double* materialProperties, int nMaterialProperties, int materialNumber )
        -> Marmot::MarmotMaterialGradientPlasticityHypoElastic* {
        return new T( materialProperties, nMaterialProperties, materialNumber );
      };
      return true;
    }

    using MaterialFactoryMap = std::unordered_map< std::string, materialFactoryFunction >;
    static MaterialFactoryMap& materialFactoryFunctionByName();
  };

} // namespace MarmotLibrary
