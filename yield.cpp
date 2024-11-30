#include <cmath>

#include "yield.hpp"

xs::polynomialYield::polynomialYield( const aceData& acefile, const uint32_t loc ) {
  const auto ncoeff = static_cast<int32_t>( acefile.xss( loc ) );
  _c.resize( ncoeff );
  for ( auto i = 0 ; i < ncoeff ; ++i ) {
    _c[i] = acefile.xss( loc+i+1 );
  }
}

double xs::polynomialYield::yield( const double erg ) const {
  auto Y = 0.0;
  for ( auto k = 0 ; k < _c.size() ; ++k ) {
    Y += _c[k] * std::pow( erg, k );
  }
  return Y;
}
