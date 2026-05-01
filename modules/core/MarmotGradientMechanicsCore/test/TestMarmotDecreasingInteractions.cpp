#include "Marmot/MarmotDecreasingInteractions.h"
#include "Marmot/MarmotTesting.h"
#include <cmath>

using namespace Marmot::Testing;
using namespace Marmot::GradientDamage::DecreasingInteractions;

void testPohTemplate()
{
  // poh(omega, eta, R) = ((1-R)*exp(-eta*omega) + R - exp(-eta)) / (1 - exp(-eta))
  // At omega=0 (no damage): poh should equal 1 for any eta and R
  throwExceptionOnFailure( checkIfEqual( poh( 0.0, 1.0, 0.0 ), 1.0 ),
                            MakeString() << __PRETTY_FUNCTION__ << ": poh(0, 1, 0) != 1" );
  throwExceptionOnFailure( checkIfEqual( poh( 0.0, 2.0, 0.3 ), 1.0 ),
                            MakeString() << __PRETTY_FUNCTION__ << ": poh(0, 2, 0.3) != 1" );

  // At omega=1, eta=1, R=0: poh = (exp(-1) - exp(-1)) / (1 - exp(-1)) = 0
  throwExceptionOnFailure( checkIfEqual( poh( 1.0, 1.0, 0.0 ), 0.0 ),
                            MakeString() << __PRETTY_FUNCTION__ << ": poh(1, 1, 0) != 0" );

  // At omega=0.5, eta=1, R=0: expected value 0.37754066...
  const double expected = ( std::exp( -0.5 ) - std::exp( -1.0 ) ) / ( 1.0 - std::exp( -1.0 ) );
  throwExceptionOnFailure( checkIfEqual( poh( 0.5, 1.0, 0.0 ), expected ),
                            MakeString() << __PRETTY_FUNCTION__ << ": poh(0.5, 1, 0) mismatch" );

  // poh must be bounded in [R, 1]: poh(omega, eta, R) >= R for all omega in [0,1]
  for ( const double omega : { 0.1, 0.3, 0.5, 0.7, 0.9 } ) {
    const double val = poh( omega, 1.0, 0.1 );
    throwExceptionOnFailure( val >= 0.1 - 1e-14,
                              MakeString() << __PRETTY_FUNCTION__ << ": poh below R at omega=" << omega );
    throwExceptionOnFailure( val <= 1.0 + 1e-14,
                              MakeString() << __PRETTY_FUNCTION__ << ": poh above 1 at omega=" << omega );
  }
}

void testFirstOrderDerivedPoh()
{
  // FirstOrderDerived::poh returns {value, dvalue/domega}
  // For omega <= 0: returns {1, 0}
  {
    auto [g, dg] = FirstOrderDerived::poh( 0.0, 1.0, 0.0 );
    throwExceptionOnFailure( checkIfEqual( g, 1.0 ),
                              MakeString() << __PRETTY_FUNCTION__ << ": value at omega=0 must be 1" );
    throwExceptionOnFailure( checkIfEqual( dg, 0.0 ),
                              MakeString() << __PRETTY_FUNCTION__ << ": derivative at omega=0 must be 0" );
  }

  // For omega > 0: value must match template function
  // dg/domega = -(1-R)*eta*exp(-eta*omega) / (1 - exp(-eta))
  {
    const double omega = 0.5;
    const double eta   = 1.0;
    const double R     = 0.0;
    auto [g, dg]       = FirstOrderDerived::poh( omega, eta, R );

    const double expected_g  = poh( omega, eta, R );
    const double expected_dg = -( 1.0 - R ) * eta * std::exp( -eta * omega ) / ( 1.0 - std::exp( -eta ) );
    throwExceptionOnFailure( checkIfEqual( g, expected_g ),
                              MakeString() << __PRETTY_FUNCTION__ << ": value mismatch at omega=0.5" );
    throwExceptionOnFailure( checkIfEqual( dg, expected_dg ),
                              MakeString() << __PRETTY_FUNCTION__ << ": derivative mismatch at omega=0.5" );
  }

  // Numerical derivative verification: dg ≈ (g(omega+h) - g(omega-h)) / (2h)
  {
    const double h     = 1e-6;
    const double omega = 0.5;
    const double eta   = 1.0;
    const double R     = 0.0;
    auto [g, dg]       = FirstOrderDerived::poh( omega, eta, R );
    const double dg_fd = ( poh( omega + h, eta, R ) - poh( omega - h, eta, R ) ) / ( 2.0 * h );
    throwExceptionOnFailure( checkIfEqual( dg, dg_fd, 1e-6 ),
                              MakeString() << __PRETTY_FUNCTION__ << ": derivative FD check failed at omega=0.5" );
  }
}

void testSecondOrderDerivedPoh()
{
  // SecondOrderDerived::poh returns {value, dvalue/domega, d2value/domega2}
  // d2g/domega2 = (1-R)*eta^2*exp(-eta*omega) / (1 - exp(-eta))
  // Note: -d2g/domega2 = -eta * dg/domega  (i.e. d2g = -eta * dg_domega)

  // At omega=0.5, eta=1, R=0: all three analytically known
  {
    const double omega  = 0.5;
    const double eta    = 1.0;
    const double R      = 0.0;
    auto [g, dg, d2g]   = SecondOrderDerived::poh( omega, eta, R );

    const double expected_g   = poh( omega, eta, R );
    const double expected_dg  = -( 1.0 - R ) * eta * std::exp( -eta * omega ) / ( 1.0 - std::exp( -eta ) );
    const double expected_d2g = ( 1.0 - R ) * eta * eta * std::exp( -eta * omega ) / ( 1.0 - std::exp( -eta ) );

    throwExceptionOnFailure( checkIfEqual( g, expected_g ),
                              MakeString() << __PRETTY_FUNCTION__ << ": value mismatch at omega=0.5" );
    throwExceptionOnFailure( checkIfEqual( dg, expected_dg ),
                              MakeString() << __PRETTY_FUNCTION__ << ": first derivative mismatch at omega=0.5" );
    throwExceptionOnFailure( checkIfEqual( d2g, expected_d2g ),
                              MakeString() << __PRETTY_FUNCTION__ << ": second derivative mismatch at omega=0.5" );
  }

  // At omega=0.5, eta=2, R=0.2
  {
    const double omega  = 0.5;
    const double eta    = 2.0;
    const double R      = 0.2;
    auto [g, dg, d2g]   = SecondOrderDerived::poh( omega, eta, R );

    const double expected_g   = poh( omega, eta, R );
    const double expected_dg  = -( 1.0 - R ) * eta * std::exp( -eta * omega ) / ( 1.0 - std::exp( -eta ) );
    const double expected_d2g = ( 1.0 - R ) * eta * eta * std::exp( -eta * omega ) / ( 1.0 - std::exp( -eta ) );

    throwExceptionOnFailure( checkIfEqual( g, expected_g ),
                              MakeString() << __PRETTY_FUNCTION__ << ": value mismatch at omega=0.5, eta=2, R=0.2" );
    throwExceptionOnFailure( checkIfEqual( dg, expected_dg ),
                              MakeString() << __PRETTY_FUNCTION__
                                           << ": first derivative mismatch at omega=0.5, eta=2, R=0.2" );
    throwExceptionOnFailure( checkIfEqual( d2g, expected_d2g ),
                              MakeString() << __PRETTY_FUNCTION__
                                           << ": second derivative mismatch at omega=0.5, eta=2, R=0.2" );
  }

  // Numerical second derivative verification: d2g ≈ (dg(omega+h) - dg(omega-h)) / (2h)
  {
    const double h     = 1e-5;
    const double omega = 0.5;
    const double eta   = 1.0;
    const double R     = 0.0;
    auto [g, dg, d2g]  = SecondOrderDerived::poh( omega, eta, R );
    auto [g_p, dg_p, d2g_p] = SecondOrderDerived::poh( omega + h, eta, R );
    auto [g_m, dg_m, d2g_m] = SecondOrderDerived::poh( omega - h, eta, R );
    const double d2g_fd = ( dg_p - dg_m ) / ( 2.0 * h );
    throwExceptionOnFailure( checkIfEqual( d2g, d2g_fd, 1e-6 ),
                              MakeString() << __PRETTY_FUNCTION__ << ": second derivative FD check failed" );
  }
}

void testFirstOrderDerivedSigmoid()
{
  // FirstOrderDerived::sigmoid returns {value, dvalue/domega}
  // For omega <= 0: returns {1, 0}
  {
    auto [g, dg] = FirstOrderDerived::sigmoid( 0.0, 0.1, 0.5, 2.0 );
    throwExceptionOnFailure( checkIfEqual( g, 1.0 ),
                              MakeString() << __PRETTY_FUNCTION__ << ": value at omega=0 must be 1" );
    throwExceptionOnFailure( checkIfEqual( dg, 0.0 ),
                              MakeString() << __PRETTY_FUNCTION__ << ": derivative at omega=0 must be 0" );
  }

  // For omega > 0: test against known analytical value
  // sigmoid(0.5, 0.1, 0.5, 2):
  //   aux1 = log(2)/log(0.5) = -1
  //   aux2 = (0.5^(-1) - 1)^2 = (2-1)^2 = 1
  //   g = 1 - (1-0.1)/(1+1) = 1 - 0.45 = 0.55
  {
    const double omega = 0.5;
    const double R     = 0.1;
    const double a     = 0.5;
    const double b     = 2.0;
    auto [g, dg]       = FirstOrderDerived::sigmoid( omega, R, a, b );
    throwExceptionOnFailure( checkIfEqual( g, 0.55, 1e-14 ),
                              MakeString() << __PRETTY_FUNCTION__ << ": value mismatch at omega=0.5" );

    // Numerical derivative verification using autodiff result at nearby points
    auto [g_p, dg_p] = FirstOrderDerived::sigmoid( omega + 1e-5, R, a, b );
    auto [g_m, dg_m] = FirstOrderDerived::sigmoid( omega - 1e-5, R, a, b );
    const double dg_fd = ( g_p - g_m ) / ( 2.0 * 1e-5 );
    throwExceptionOnFailure( checkIfEqual( dg, dg_fd, 1e-5 ),
                              MakeString() << __PRETTY_FUNCTION__ << ": derivative FD check failed at omega=0.5" );
  }

  // At omega=1, R=0.1, a=0.5, b=2: sigmoid(1,0.1,0.5,2) = R = 0.1
  // (at omega=1, aux1=-1, pow(1,-1)-1=0, aux2=0^2=0, g=1-0.9/1=0.1)
  {
    auto [g, dg] = FirstOrderDerived::sigmoid( 1.0, 0.1, 0.5, 2.0 );
    throwExceptionOnFailure( checkIfEqual( g, 0.1 ),
                              MakeString() << __PRETTY_FUNCTION__ << ": sigmoid(1,0.1,0.5,2) != R=0.1" );
  }
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ testPohTemplate,
                                                       testFirstOrderDerivedPoh,
                                                       testSecondOrderDerivedPoh,
                                                       testFirstOrderDerivedSigmoid };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
