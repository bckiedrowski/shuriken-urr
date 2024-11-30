#include <cmath>

#include "cell.hpp"

bool geom::cell::within( const util::point p ) const {
  //loop over surfaces, if the sign of the evaulation of each surface matches
  //with each sense, then we are in the surface; otherwise not
  for ( auto i = 0 ; i < _surfaces.size() ; ++i ) {
    if ( std::copysign( 1, _surfaces[i]->eval(p) ) != _senses[i] ) return false;
  }
  return true;
}

geom::nextSurface geom::cell::next_surface( const util::point p, const util::point u ) const {
  //find surface with minimum positive distance
  nextSurface next;
  for ( auto s : _surfaces ) {
    const auto dist = s->distance( p, u );
    if ( dist < next.distance ) {
      next.distance = dist;
      next.surface  = s;
    }
  }
  return next;
}

std::shared_ptr< geom::cell > geom::find_cell( const util::point p,
  const std::vector< std::shared_ptr< cell > > cells ) {
  //search a vector of cells for the first cell that the point p is within
  //if no cell found, return the null pointer
  for ( auto c : cells ) {
    if ( c->within(p) ) { return c; }
  }
  return nullptr;
}
