#include <cmath>
#include <limits>

#include "quadSolve.hpp"

#include "surface.hpp"

//---------------------------------------------------------------------------------------
util::point geom::surface::normal( const util::point p ) const {
  util::point g = grad( p );
  return g / std::sqrt( util::dot_product(g,g) );
}

util::point geom::surface::reflect( const util::point p, const util::point u ) const {
  util::point g = grad( p );
  return u - g * 2.0 * dot_product( u, g ) / dot_product( g, g );
}

//---------------------------------------------------------------------------------------
double geom::plane::eval( const util::point p ) const {
  return a*p.x + b*p.y + c*p.z - d;
}

double geom::plane::distance( const util::point p, const util::point u ) const {
  const double dist = (d - a * p.x - b * p.y - c * p.z) / (a * u.x + b * u.y + c * u.z);
  return dist >= 0.0 ? dist : std::numeric_limits<double>::max();
}

util::point geom::plane::grad( const util::point p ) const {
  return util::point( a, b, c );
}

//---------------------------------------------------------------------------------------
double geom::sphere::eval( const util::point p ) const {
//  return std::pow(p.x - x0, 2) + std::pow(p.y - y0, 2) + std::pow(p.z - x0, 2) - r*r;
  const util::point q(p.x - x0, p.y - y0, p.z - z0);
  return q.x*q.x + q.y*q.y + q.z*q.z - r*r;
}

double geom::sphere::distance( const util::point p, const util::point u ) const {
  const util::point q(p.x - x0, p.y - y0, p.z - z0);
  const double b = 2.0 * (q.x * u.x + q.y * u.y + q.z * u.z);
  const double c = eval(p);

  return util::quad_solve(1.0, b, c);
}

util::point geom::sphere::grad( const util::point p ) const {
  return util::point( 2.0*(p.x - x0), 2.0*(p.y - y0), 2.0*(p.z - z0) );
}

//---------------------------------------------------------------------------------------
double geom::cylinder_z::eval( const util::point p ) const {
  return std::pow(p.x - x0, 2) + std::pow(p.y - y0, 2) - r*r;
}

double geom::cylinder_z::distance( const util::point p, const util::point u ) const {
  const util::point q(p.x - x0, p.y - y0, p.z);

  const double a = std::pow(u.x, 2) + std::pow(u.y, 2);
  const double b = 2.0 * (q.x * u.x + q.y * u.y);
  const double c = eval(p);

  return a == 0 && b == 0 ? 0.0 : util::quad_solve(a, b, c);
}

util::point geom::cylinder_z::grad( const util::point p ) const {
  return util::point( 2.0*(p.x - x0), 2.0*(p.y - y0), 0.0 );
}
