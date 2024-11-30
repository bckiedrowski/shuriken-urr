#include <cmath>

#include "binarySearch.hpp"

#include "angularDistribution.hpp"

double xs::isotropicAngDist::sample( util::RNG::ptr rng ) const {
  return 2.0 * rng->sample() - 1.0;
}

xs::equiprobable32BinAngDist::equiprobable32BinAngDist(
  const aceData& acefile, const uint32_t loc ) {
  for ( auto i = 0 ; i < 33 ; ++i ) {
    _mu[i] = acefile.xss( loc + i );
  }
}

double xs::equiprobable32BinAngDist::sample( util::RNG::ptr rng ) const {
  const auto r = rng->sample();
  const auto k = std::floor( 32*r );
  return _mu[k] + ( 32*r - k )*( _mu[k+1] - _mu[k] );
}

xs::tabularAngDist::tabularAngDist(
  const aceData& acefile, const uint32_t loc ) {

  _interp_type = static_cast<interp_t>( acefile.xss( loc ) );
  assert( _interp_type == interp_t::HISTOGRAM || _interp_type == interp_t::LINLIN );

  _ncos = static_cast<int32_t>( acefile.xss( loc+1 ) );

  _mu.resize(  _ncos );
  _pdf.resize( _ncos );
  _cdf.resize( _ncos );
  for ( auto i = 0 ; i < _ncos ; ++i ) {
    _mu[i]  = acefile.xss( loc + 2 + i );
    _pdf[i] = acefile.xss( loc + 2 +   _ncos + i );
    _cdf[i] = acefile.xss( loc + 2 + 2*_ncos + i );
  }
}

double xs::tabularAngDist::sample( util::RNG::ptr rng ) const {
  //search cdf grid using random number
  //the cdf grid is often small, so use a linear search when this is the case
  //and a binary search when it is not
  const auto r = rng->sample();
  uint32_t   k = 0;
  if ( _cdf.size() <= 8 ) {
    for ( ; k < _cdf.size()-1 ; ++k ) {
      if ( r < _cdf[k+1] ) break;
    }
  }
  else {
    k = util::binary_search( r, _cdf );
  }

  //return value based on interpolation type
  if ( _interp_type == interp_t::HISTOGRAM ) {
    return _mu[k] + ( r - _cdf[k] )/_pdf[k];
  }
  else { //lin-lin interpolation
    const auto f = ( _pdf[k+1] - _pdf[k] )/( _mu[k+1] - _mu[k] );
    if ( std::fabs(f) > 1.0e-9 ) {
      return _mu[k] + ( std::sqrt( _pdf[k]*_pdf[k] + 2*f*(r - _cdf[k]) ) - _pdf[k] )/f;
    }
    else {
      //pdf is flat (or effectively so) in this range and the formula above is invalid
      //and a histogram treatment is correct or a good approximation
      return _mu[k] + ( r - _cdf[k] )/_pdf[k];
    }
  }
}

xs::angularDistributionTable::angularDistributionTable( const aceData& acefile, const uint32_t mt_offset ) {

  const auto loc_land   = acefile.jxs(8);  //location of angular distribution locators
  const auto loc_and    = acefile.jxs(9);  //location of angular distribution block

  const auto locb = static_cast<int32_t>( acefile.xss( loc_land + mt_offset ) );
  if ( locb == 0 ) {
    //isotropic scattering over entire energy range
    _present   = true;
    _isotropic = true;
  }
  else if ( locb == -1 ) {
    //correlated angular distribution given elsewhere (via ACE law 44 or 61)
  }
  else {
    //angular distributions provided as function of energy
    _present = true;
    const auto loc_adist = loc_and + locb - 1;

    _nerg  = static_cast<int32_t>( acefile.xss( loc_adist ) );
    _egrid.resize( _nerg );
    _ang.resize( _nerg );
    for ( auto i = 0 ; i < _nerg ; ++i ) {
      _egrid[i] = acefile.xss( loc_adist + i + 1 );
      const auto lc = static_cast<int32_t>( acefile.xss( loc_adist + _nerg + i + 1 ) );
      if ( lc == 0 ) {
        //isotropic
        _ang[i] = std::make_shared< isotropicAngDist >();
      }
      else if ( lc > 0 ) {
        //32 equiprobable bins
        _ang[i] = std::make_shared< equiprobable32BinAngDist >( acefile, loc_and + lc - 1 );
      }
      else {
        //tabular
        _ang[i] = std::make_shared< tabularAngDist >( acefile, loc_and + std::abs(lc) - 1 );
      }
    }
  }
}

double xs::angularDistributionTable::sample( const double erg, util::RNG::ptr rng ) const {

  assert( _present ); //should not call sample routine if table not present (host function needs to test)

  if ( _isotropic ) {
    return 2.0 * rng->sample() - 1.0;
  }
  else {
    //search energy grid and randomly sample the table
    //check out of bounds, sample table at endpoint
    if ( erg <  _egrid[0] )       { return _ang[0]->sample( rng ); }
    if ( erg >= _egrid[_nerg-1] ) { return _ang[_nerg-1]->sample( rng ); }

    //binary search on grid of values
    const auto k = util::binary_search( erg, _egrid );

    //randomly sample table based on proximity in energy to the table (lin-lin interpolation)
    const auto r = ( erg - _egrid[k] )/( _egrid[k+1] - _egrid[k] );
    return rng->sample() >= r ? _ang[k]->sample( rng ) : _ang[k+1]->sample( rng );
  }

}
