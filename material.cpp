#include "material.hpp"

void mc::material::add( const double aden, std::shared_ptr<xs::neutronData> xs  ) {
  _aden.push_back( aden );
  _xs.push_back( xs );
  _id.push_back( xs->id() );
  _size++;
}

void mc::material::fill_total_xs( const double erg, xs::xsBuffer& xs_buffer, util::RNG::ptr rng ) const {
  xs_buffer.total_macro_xs = 0.0;
  for ( auto i = 0 ; i < _size ; ++i ) {
    const auto k = _xs[i]->id();
    if ( xs_buffer.erg[k] != erg ) {
      _xs[i]->fill_buffer( erg, xs_buffer, rng );
    }
    xs_buffer.total_macro_xs += _aden[i] * xs_buffer.total_xs[k];
  }
}

void mc::material::fill_nufission_xs( const double erg, xs::xsBuffer& xs_buffer, util::RNG::ptr rng ) const {
  xs_buffer.nufission_macro_xs = 0.0;
  for ( auto i = 0 ; i < _size ; ++i ) {
    const auto k = _xs[i]->id();
    if ( xs_buffer.erg[k] != erg ) {
      _xs[i]->fill_buffer( erg, xs_buffer, rng );
    }
    xs_buffer.nufission_macro_xs += _aden[i] * xs_buffer.nufission_xs[k];
  }
}

//sample an isotope undergoing a specific reaction (usually collision or fission)
//for speed need to pass in the sum weighted by the atomic density of each isotope
//and a vector of partial microscopic cross sections of those isotopes (normally from xs_buffer)
//note: the vector of partials must be ordered by global indexing which usually differs
//      from order within material definition
std::shared_ptr< xs::neutronData > mc::material::select_isotope( const double sum_macro,
  const std::vector<double> partial_micro, util::RNG::ptr rng ) const {
  double  t = rng->sample() * sum_macro;
  int32_t i = -1;
  while ( t > 0 ) {
    ++i;
    const auto k = _id[i];
    t -= _aden[i] * partial_micro[k];
  }
  return _xs[i];
}
