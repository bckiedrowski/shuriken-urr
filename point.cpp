
#include "point.hpp"

util::point util::point::operator+( const point& p ) const {
  return point( x + p.x, y + p.y, z + p.z );
}

util::point util::point::operator-( const point& p ) const {
  return point( x - p.x, y - p.y, z - p.z );
}

util::point util::point::operator*( const double a ) const {
  return point( a*x, a*y, a*z );
}

util::point util::point::operator/( const double a ) const {
  return point( x/a, y/a, z/a );
}

double util::dot_product( const util::point& p, const util::point& q ) {
  return p.x*q.x + p.y*q.y + p.z*q.z;
}
