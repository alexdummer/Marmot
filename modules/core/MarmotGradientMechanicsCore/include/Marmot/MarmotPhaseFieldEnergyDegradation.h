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

#include <tuple>

namespace Marmot::PhaseField {

  namespace EnergyDegradationFunctions {
    // degradation functions acc. Kuhn, Schlüter, Müller (2015)
    template < typename T >
    T quadratic( const T pf )
    {
      return pow( 1. - pf, 2 );
    }

    template < typename T >
    T cubic( const T pf )
    {
      const T s = 1. - pf;
      return 3. * pow( s, 2 ) + 2. * pow( s, 3 );
    }

    template < typename T >
    T quartic( const T pf )
    {
      const T s = 1. - pf;
      return 4 * pow( s, 3 ) + 3 * pow( s, 4 );
    }

    // generic degradation function acc. Wu, JMPS (2017)
    template < typename T >
    T generic( const T pf, const double p, const double a1, const double a2, const double a3 )
    {
      const T s      = 1. - pf;
      const T P      = 1. + a2 * pf * ( 1. + a3 * pf );
      const T Q      = a1 * pf * P;
      const T result = pow( s, p ) / ( pow( s, p ) + Q );

      return result;
    }

    namespace SecondOrderDerived {
      std::tuple< double, double, double > quadratic( const double pf );
      std::tuple< double, double, double > cubic( const double pf );
      std::tuple< double, double, double > quartic( const double pf );
      std::tuple< double, double, double > generic( const double pf,
                                                    const double p,
                                                    const double a1,
                                                    const double a2,
                                                    const double a3 );
    }; // namespace SecondOrderDerived
  };   // namespace EnergyDegradationFunctions
};     // namespace Marmot::PhaseField
