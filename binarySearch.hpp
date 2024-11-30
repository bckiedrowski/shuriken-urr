#ifndef _UTIL_BINARY_SEARCH_HEADER_
#define _UTIL_BINARY_SEARCH_HEADER_

#include <cassert>
#include <cstdint>
#include <vector>

namespace util {
  //binary search returns lower index
  uint32_t binary_search( const double x, const std::vector<double>& grid, const uint32_t start = 0 );
  uint32_t binary_search( const double x, const std::vector<double>& grid, const uint32_t start, const uint32_t end );
}

#endif
