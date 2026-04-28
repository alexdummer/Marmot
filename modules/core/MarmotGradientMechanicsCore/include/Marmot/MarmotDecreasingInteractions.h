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

    template < typename T >
    T sigmoid( T omega, double R, double a, double b )
    {
      /* auto    sig = [&]( const T om ) { return 1. / ( 1. + exp( -eta * ( om - a ) ) ); }; */
      /* const T s   = ( 1. - R ) / ( sig( 1. ) - sig( 0. ) ); */

      /* return 1. + s * ( sig( 0 ) - sig( omega ) ); */

      const double aux1 = log( 2 ) / log( a );
      const T      aux2 = pow( eval( pow( omega, aux1 ) - 1. ), b );

      return 1. - ( 1. - R ) / ( 1. + aux2 );
    }
    namespace FirstOrderDerived {

      std::tuple< double, double > poh( double omega, double eta, double R );

      std::tuple< double, double > sigmoid( double omega, double R, double a, double b );

    } // namespace FirstOrderDerived

    namespace SecondOrderDerived {

      std::tuple< double, double, double > poh( double omega, double eta, double R );

      std::tuple< double, double, double > sigmoid( double omega, double R, double a, double b );

    } // namespace SecondOrderDerived
  }   // namespace DecreasingInteractions

} // namespace Marmot::GradientDamage
