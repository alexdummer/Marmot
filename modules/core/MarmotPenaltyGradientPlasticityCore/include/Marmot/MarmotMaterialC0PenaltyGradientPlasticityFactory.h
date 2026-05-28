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
#pragma once
#include "Marmot/MarmotMaterialC0PenaltyGradientPlasticity.h"
#include <cassert>
#include <functional>
#include <string>
#include <unordered_map>

namespace MarmotLibrary {

  /**
   * @brief Factory class for creating MarmotMaterialC0PenaltyGradientPlasticity instances.
   *
   * Uses the Factory Design Pattern to manage creation of different
   * C0 penalty gradient plasticity materials based on unique material names.
   */
  class MarmotMaterialC0PenaltyGradientPlasticityFactory {
  public:
    using materialFactoryFunction = std::function<
      MarmotMaterialC0PenaltyGradientPlasticity*( const double* materialProperties,
                                                  int           nMaterialProperties,
                                                  int           materialNumber ) >;

    MarmotMaterialC0PenaltyGradientPlasticityFactory() = delete;

    /**
     * @brief Create a material instance based on its name and properties.
     * @param[in] materialName Unique name of the material to create.
     * @param[in] materialProperties Array of material properties.
     * @param[in] nMaterialProperties Number of properties in the array.
     * @param[in] materialNumber Unique identifier for the material instance.
     * @return Pointer to the created material instance.
     */
    static MarmotMaterialC0PenaltyGradientPlasticity* createMaterial( const std::string& materialName,
                                                                      const double*      materialProperties,
                                                                      int                nMaterialProperties,
                                                                      int                materialNumber );

    /**
     * @brief Register a material with its factory function.
     * @tparam T The concrete material class to register.
     * @param[in] materialName Name of the material.
     * @return True if registration was successful.
     */
    template < class T >
    static bool registerMaterial( const std::string& materialName )
    {
      auto& map = materialFactoryFunctionByName();

      assert( map.find( materialName ) == map.end() && "Material already registered!" );

      map[materialName] = []( const double* materialProperties,
                              int           nMaterialProperties,
                              int           materialNumber ) -> MarmotMaterialC0PenaltyGradientPlasticity* {
        return new T( materialProperties, nMaterialProperties, materialNumber );
      };
      return true;
    }

    using MaterialFactoryMap = std::unordered_map< std::string, materialFactoryFunction >;
    /// @brief Get the map of material factory functions by material name.
    static MaterialFactoryMap& materialFactoryFunctionByName();
  };
} // namespace MarmotLibrary
