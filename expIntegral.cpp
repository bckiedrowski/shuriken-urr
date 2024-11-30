#include "expIntegral.hpp"

double util::expint( uint32_t n, double x ) {
  if ( n == 1 ) {
    return -std::expint(-x);
  }
  else {
    return ( std::exp(-x) - x * util::expint( n-1, x ) )/(n-1);
  }
}
