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
#include <stdexcept>
#include <string>

namespace Marmot {

  /** @class StressUpdateFailed
   *  @brief Exception thrown when a stress update (e.g., return-mapping iteration) fails to converge
   *         inside a material model.
   *
   *  Inherits from std::runtime_error so that existing catch(std::runtime_error&) handlers still
   *  work. Callers that need to distinguish a material-model failure from other runtime errors
   *  (e.g., I/O errors originating in third-party libraries) should catch StressUpdateFailed
   *  explicitly.
   */
  class StressUpdateFailed : public std::runtime_error {
  public:
    explicit StressUpdateFailed( const std::string& message ) : std::runtime_error( message ) {}
  };

} // namespace Marmot
