#include <cmath>

#include "binarySearch.hpp"

#include "interpGrid.hpp"

xs::interpGrid::interpGrid( const aceData& acefile, const uint32_t loc ) {
  const auto num_interp = static_cast<int32_t>( acefile.xss(loc) );
  const auto num_vals   = static_cast<int32_t>( acefile.xss(loc + 2*num_interp + 1) );

  //set up coarse grid
  if ( num_interp == 0 ) {
    //special case of a single bin with LINLIN
    _interp_bounds.resize(1);
    _interp_bounds[0] = num_vals-1;

    _interp_types.resize(1);
    _interp_types[0] = interp_t::LINLIN;
  }
  else {
    _interp_bounds.resize(num_interp);
    _interp_types.resize(num_interp);
    for ( auto i = 0 ; i < num_interp ; ++i ) {
      _interp_bounds[i] = static_cast<int32_t>( acefile.xss(loc + i + 1) - 1 );
      _interp_types[i]  = static_cast<interp_t>( acefile.xss(loc + num_interp + i + 1) );
    }
  }

  const auto loc_xgrid = loc + 2*num_interp + 2;
  const auto loc_ygrid = loc_xgrid + num_vals;

  _xgrid.resize( num_vals );
  _ygrid.resize( num_vals );
  for ( auto i = 0 ; i < num_vals ; ++i ) {
    _xgrid[i] = acefile.xss( loc_xgrid + i );
    _ygrid[i] = acefile.xss( loc_ygrid + i );
  }
}

double xs::interpGrid::operator()( const double x ) const {
  //check out of bounds, return value endpoint
  if ( x <  _xgrid[0] ) { return _ygrid[0]; }
  if ( x >= _xgrid[_xgrid.size()-1] ) { return _ygrid[_ygrid.size()-1]; }

  //often the xgrid is very small and binary search is more expensive, so perform
  //a linear search if the grid is small and a binary search if it is not
  uint32_t ibot = 0;
  if ( _xgrid.size() <= 8 ) {
    for ( ; ibot < _xgrid.size()-1 ; ++ibot ) {
      if ( x < _xgrid[ibot+1] ) break;
    }
  }
  else {
    ibot = util::binary_search( x, _xgrid );
  }
  uint32_t itop = ibot + 1;

  //grab interpolation scheme from coarse grid with linear search since there
  //are usually very few (usually one) grid point
  uint32_t k = 0;
  for ( k = 0 ; k < _interp_bounds.size() ; ++k ) {
    if ( itop <= _interp_bounds[k] ) break;
  }
  const auto interp = _interp_types[k];

  //perform interpolation
  auto r = 0.0;
  switch( interp ) {
  case interp_t::HISTOGRAM:
    //return value at lower limit as specified by ENDF-6 format
    return _ygrid[ibot];
  case interp_t::LINLIN:
    r = ( x - _xgrid[ibot] )/( _xgrid[itop] - _xgrid[ibot] );
    return _ygrid[ibot] + r * ( _ygrid[itop] - _ygrid[ibot] );
  case interp_t::LINLOG:
    r = std::log(x/_xgrid[ibot]) / std::log(_xgrid[itop]/_xgrid[ibot]);
    return _ygrid[ibot] + r * ( _ygrid[itop] - _ygrid[ibot] );
  case interp_t::LOGLIN:
    r = ( x - _xgrid[ibot] )/( _xgrid[itop] - _xgrid[ibot] );
    return _ygrid[ibot] * std::pow( _ygrid[itop]/_ygrid[ibot], r );
  case interp_t::LOGLOG:
    r = std::log(x/_xgrid[ibot]) / std::log(_xgrid[itop]/_xgrid[ibot]);
    return _ygrid[ibot] * std::pow( _ygrid[itop]/_ygrid[ibot], r );
  default:
    throw("unknown interpolation type");
    return 0.0;
  }
}
