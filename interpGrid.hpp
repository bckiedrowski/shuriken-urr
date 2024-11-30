#ifndef _XS_INTERP_GRID_HEADER_
#define _XS_INTERP_GRID_HEADER_

#include <vector>

#include "ace.hpp"

namespace xs {

  enum class interp_t {
    HISTOGRAM = 1, LINLIN = 2, LINLOG = 3, LOGLIN = 4, LOGLOG = 5 };
  // ex: LINLOG is read as y is linear in log of x

  class interpGrid {
    private:
      //1d data on a fine grid
      std::vector<double>   _xgrid;
      std::vector<double>   _ygrid;

      //interpolation schemes on a coarse grid containing upper bound indices
      std::vector<uint32_t> _interp_bounds;
      std::vector<interp_t> _interp_types;
    public:
      //construct from section of ace file
      interpGrid( const aceData& acefile, const uint32_t loc );

      double operator()( const double x ) const;
  };

}

#endif
