
#include "binarySearch.hpp"

//version with optional third start argument
uint32_t util::binary_search( const double x, const std::vector<double>& grid, const uint32_t start ) {
  //assumes x is within range, exclusive of top grid point
  assert( x >= grid[start] && x < grid[grid.size()-1] );

  uint32_t ibot = start;
  uint32_t itop = grid.size()-1;
  while ( itop - ibot > 1 ) {
    uint32_t imid = (ibot + itop)/2;
    if ( x < grid[imid] ){
      itop = imid;
    }
    else {
      ibot = imid;
    }
  }
  return ibot;
}

//version where both third start and bottom end arguments MUST be specified
//this instance is normally called after use of the hash functions
uint32_t util::binary_search( const double x, const std::vector<double>& grid, const uint32_t start, const uint32_t end ) {
  //assumes x is within range, exclusive of top grid point
  assert( x >= grid[start] && x < grid[end] );

  uint32_t ibot = start;
  uint32_t itop = end;
  while ( itop - ibot > 1 ) {
    uint32_t imid = (ibot + itop)/2;
    if ( x < grid[imid] ){
      itop = imid;
    }
    else {
      ibot = imid;
    }
  }
  return ibot;
}
