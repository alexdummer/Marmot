#include "Marmot/MarmotDecreasingInteractions.h"
#include "autodiff/forward/dual.hpp"

namespace Marmot::GradientDamage::DecreasingInteractions {

  namespace FirstOrderDerived {

    std::tuple< double, double > poh( double omega, double eta, double R )
    {

      const double g         = Marmot::GradientDamage::DecreasingInteractions::poh( omega, eta, R );
      const double dg_domega = -( 1. - R ) * eta * exp( -eta * omega ) / ( 1. - exp( -eta ) );

      /* double dg_domega = 0; */
      /* double g         = 1; */
      /* if ( omega > 0 ) { */
      /*   autodiff::dual omega_dual( omega ); */
      /*   const auto     g_dual = Marmot::GradientDamage::DecreasingInteractions::poh( omega_dual, eta, R ); */
      /*   g                     = g_dual.val; */
      /*   dg_domega             = g_dual.grad; */
      /* } */
      if ( omega > 0 ) {
        return { g, dg_domega };
      }
      else {
        return { 1, 0 };
      }
    }
    std::tuple< double, double > sigmoid( double omega, double R, double a, double b )
    {
      /* auto         sig = [&]( const double om ) { return 1. / ( 1. + exp( -eta * ( om - a ) ) ); }; */
      /* const double s   = ( 1. - R ) / ( sig( 1. ) - sig( 0. ) ); */

      /* const double dg_domega =  - s * eta * exp( -eta * ( omega - a ) ) / ( 1. + exp( -eta * ( omega - a ) ) )  /
       * ( 1. + exp( -eta * ( omega - a ) ) ); */
      double dg_domega = 0;
      double g         = 1;
      if ( omega > 0 ) {
        autodiff::dual omega_dual( omega );
        const auto     g_dual = Marmot::GradientDamage::DecreasingInteractions::sigmoid( omega_dual, R, a, b );
        g                     = g_dual.val;
        dg_domega             = g_dual.grad;
        /* const double aux1 = log( 2 ) / log( a ); */
        /* const double aux2 = pow( pow( omega, aux1 ) - 1., b ); */
        /* g =  1. - (1. - R ) / ( 1. + aux2 ); */
        /* dg_domega = ( 1 - R ) * pow( 1 + aux2, -2 ) * b * pow( pow( omega, aux1 ) - 1 ,  b - 1 ) * aux1 * pow( omega,
         * aux1 - 1 ) ; */
      }
      return { g, dg_domega };
    }
  } // namespace FirstOrderDerived
  namespace SecondOrderDerived {

    std::tuple< double, double, double > poh( double omega, double eta, double R )
    {
      const double g           = Marmot::GradientDamage::DecreasingInteractions::poh( omega, eta, R );
      const double dg_domega   = -( 1. - R ) * eta * exp( -eta * omega ) / ( 1. - exp( -eta ) );
      const double d2g_domega2 = ( 1. - R ) * eta * eta * exp( -eta * omega ) / ( 1. - exp( -eta ) );
      return { g, dg_domega, d2g_domega2 };
    }

    std::tuple< double, double, double > sigmoid( double omega, double R, double a, double b )
    {
      throw std::runtime_error( "Second order derivative of sigmoid not implemented" );
    }
  } // namespace SecondOrderDerived

} // namespace Marmot::GradientDamage::DecreasingInteractions
