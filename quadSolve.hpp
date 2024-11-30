#ifndef _UTIL_QUAD_SOLVE_HEADER_
#define _UTIL_QUAD_SOLVE_HEADER_

#include <vector>

namespace util {
  //quadratic equation solver returns smallest real positive root
  //if one does not exist, returns largest possible double
  double quad_solve( const double a, const double b, const double c );
}

#endif
