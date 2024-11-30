#include <cmath>
#include <iomanip>
#include <iostream>

#include "binarySearch.hpp"

#include "estimator.hpp"

mc::estimator::estimator( const std::string name, const int32_t id, const mc::estimator_type_t type,
                          const std::optional< std::set<uint32_t> >  geom_id,
                          const std::optional< uint32_t >            zaid,
                          const std::optional< xs::endf_mt_t >       mt,
                          const std::optional< std::vector<double> >& ebins ) {
  _name = name;
  _id   = id;
  _type = type;
  if ( geom_id ) _geom_id = geom_id.value();
  if ( zaid    ) _zaid    = zaid.value();
  if ( mt      ) _mt      = mt.value();
  if ( ebins   ) {
    _egrid = ebins.value();
    assert( _egrid.size() >= 2 );
    _nebins = _egrid.size()-1;
  }

  if ( _mt == xs::endf_mt_t::TOTAL     || _mt == xs::endf_mt_t::ELASTIC ||
       _mt == xs::endf_mt_t::FISSION   || _mt == xs::endf_mt_t::Z_GAMMA ||
       _mt == xs::endf_mt_t::NUFISSION ) {
    _within_xs_buffer = true;
  }

  //element zero is the sum total over all the energy bins
  //if energy bins not specified, estimator is over all energies, allocate one element
  //if a single energy bin is specified, only allocate one element
  //if multiplie energy bins specified, allocate a space for the total and one for each bin
  const auto array_sizes = _nebins > 1 ? 1 + _nebins : 1;
  assert( array_sizes != 2 );
  _sum_hist.resize( array_sizes, 0.0 );
  _mean.resize( array_sizes, 0.0 );
  _stdv.resize( array_sizes, 0.0 );
}

bool mc::estimator::within_geom_set( const uint32_t geom_id ) const {
  if ( ! _geom_id.empty() ) {
    return _geom_id.find(geom_id) != _geom_id.end();
  }
  return true;
}

bool mc::estimator::within_energy_range( const double erg ) const {
  if ( _nebins > 0 ) {
    if ( erg < _egrid[0] || erg >= _egrid[_nebins] ) return false;
  }
  return true;
}

uint32_t mc::estimator::energy_bin( const double erg ) const {
  assert( within_energy_range( erg ) );
  if ( _nebins >= 2 ) {
    if ( _nebins < 10 ) {
      //linear search for small energy grid
      uint32_t ie = 0;
      for ( ie = 1 ; ie < _nebins ; ++ie ) {
        if ( _egrid[ie] > erg ) break;
      }
      return ie;
    }
    else {
      return 1 + util::binary_search( erg, _egrid );
    }
  }
  return 0;   //return 0 for a single bin
}

void mc::estimator::score( const double contribution, const uint32_t geom_id, const double erg,
                           const std::shared_ptr< geom::cell > cell_ptr,
                           const xs::xsBuffer& xs_buffer ) {
  //filter out non-matching geometry IDs
  if ( ! _geom_id.empty() && _geom_id.find(geom_id) == _geom_id.end() ) return;

  //filter out if not within energy range and grab index is within
  //note that zero index implies a single energy bin
  if ( ! within_energy_range( erg ) ) return;
  const auto ie = energy_bin( erg );

  //if no specific reaction type, just score the base response and leave
  if ( _mt == xs::endf_mt_t::NONE ) {
    _sum_hist[0]  += contribution;
    if ( ie > 0 ) _sum_hist[ie] += contribution;
    return;
  }

  //otherwise need to have cross section multiplier
  const auto mat = cell_ptr->material();
  if ( _zaid == 0 ) {
    //all isotopes in material
    for ( auto i = 0 ; i < mat->size() ; ++i ) {
      const auto iso = mat->xs(i);
      const auto k   = iso->id();
      const auto s   = contribution * mat->aden(i) *
        ( _within_xs_buffer ? xs_buffer( _mt, k ) : iso->find_xs( _mt, erg ) );
      _sum_hist[0] += s;
      if ( ie > 0 ) _sum_hist[ie] += s;
    }
  }
  else {
    //specific isotope
    std::shared_ptr< xs::neutronData > iso;
    uint32_t i;
    bool     found = false;
    for ( i = 0 ; i < mat->size() ; ++i ) {
      iso = mat->xs(i);
      if ( iso->zaid() == _zaid ) { found = true; break; }
    }
    if ( found ) {
      const auto k  = iso->id();
      const auto s   = contribution * mat->aden(i) *
        ( _within_xs_buffer ? xs_buffer( _mt, k ) : iso->find_xs( _mt, erg ) );
      _sum_hist[0] += s;
      if ( ie > 0 ) _sum_hist[ie] += s;
    }
  }
}

void mc::estimator::end_history() {
  for ( auto i = 0 ; i < _sum_hist.size() ; ++i ) {
    _mean[i] += _sum_hist[i];
    _stdv[i] += _sum_hist[i] * _sum_hist[i];
    _sum_hist[i] = 0.0;
  }
}

double mc::estimator::mean( const double normalization, const uint32_t k ) const {
  return _mean[k] / normalization;
}

double mc::estimator::stdv( const double normalization, const uint32_t k ) const {
  const auto m1 = mean( normalization, k );
  return std::sqrt( (_stdv[k]/normalization - m1*m1 )/normalization );
}

void mc::estimator::write( const double normalization ) const {
  if ( _nebins == 0 ) {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << _name << ' ' << std::setw(12) << mean( normalization, 0 )
              << " +/- " << stdv( normalization, 0 ) << '\n';
  }
  else {
    //write each energy bin
    for ( auto k = 1 ; k < _mean.size() ; ++k ) {
      std::cout << _name << ' ' << std::scientific << std::setprecision(4) << std::setw(10)
                << _egrid[k-1] << " to " << _egrid[k] << " = "
                << std::setprecision(6) << std::setw(12)
                << mean( normalization, k ) << " +/- " << stdv( normalization, k ) << '\n';
    }
    //write total
    std::cout << _name << ' ' << std::scientific << std::setprecision(6);
    std::cout << " Total = " << std::setw(12) << mean( normalization, 0 ) << " +/- " << stdv( normalization, 0 ) << '\n';
  }
  std::cout << std::fixed;
}

void mc::score_estimators_of_type( const estimator_type_t type,
    const std::vector< std::shared_ptr<estimator> >& estimators,
    const double contribution, const uint32_t geom_id, const double erg,
    const std::shared_ptr< geom::cell > cell_ptr, const xs::xsBuffer& xs_buffer ) {

  for ( auto est : estimators ) {
    if ( est->type() == type ) est->score( contribution, geom_id, erg, cell_ptr, xs_buffer );
  }

}

void mc::end_history_estimators( const std::vector< std::shared_ptr<estimator> >& estimators ) {
  for ( auto est : estimators ) {
    est->end_history();
  }
}

void mc::write_estimators( const std::vector< std::shared_ptr<estimator> >& estimators, const double normalization ) {
  for ( auto est : estimators ) {
    est->write( normalization );
  }
/*
  size_t max_width = 0;
  for ( auto est : estimators ) {
    max_width = std::max( max_width, est->name().length() );
  }
  for ( auto est : estimators ) {
    const auto name = est->name();
    const auto mean = est->mean( normalization );
    const auto stdv = est->stdv( normalization );
    std::cout << std::scientific << std::setprecision(6);
    std::cout << std::setw(max_width+2) << name << ' ' << std::setw(12) << mean << " +/- " << stdv << '\n';
  }
  std::cout << std::fixed;
*/
}
