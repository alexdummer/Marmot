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
 * Alexander Dummer alexander.dummer@uibk.ac.at
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
#include <cmath>
#include <tuple>

namespace Marmot::GradientDamage {

  namespace DecreasingInteractions {

    template < typename T >
    T poh( T omega, double eta, double R )
    {
      const T g = ( ( 1. - R ) * exp( -eta * omega ) + R - exp( -eta ) ) / ( 1. - exp( -eta ) );
      return g;
    }

    namespace FirstOrderDerived {

      std::tuple< double, double > poh( double omega, double eta, double R );

    } // namespace FirstOrderDerived

    namespace SecondOrderDerived {

      std::tuple< double, double, double > poh( double omega, double eta, double R );

    } // namespace SecondOrderDerived
  }   // namespace DecreasingInteractions

} // namespace Marmot::GradientDamage
