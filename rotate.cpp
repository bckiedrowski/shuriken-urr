#include <cmath>

#include "rotate.hpp"

//perform coordinate rotation on dir by direction cosine mu and azimuthal angle azi
util::point util::rotate( const util::point& dir, const double mu, const double azi ) {
  const auto sin_azi = std::sin(azi);
  const auto cos_azi = std::cos(azi);
  const auto a = std::sqrt( std::fmax( 0.0, 1.0 - mu*mu ) );
  auto b = std::sqrt( std::fmax( 0.0, 1.0 - dir.z*dir.z ) );

  //check divide by zero and perform calculation in along another axis if so
  double x, y, z;
  if ( b > 1.e-10 ) {
    x = mu*dir.x + a*(dir.x * dir.z * cos_azi - dir.y * sin_azi)/b;
    y = mu*dir.y + a*(dir.y * dir.z * cos_azi + dir.x * sin_azi)/b;
    z = mu*dir.z - a*b*cos_azi;
  }
  else {
    b  = std::sqrt(1.0 - dir.y * dir.y );
    x = mu*dir.x + a*(dir.x * dir.y * cos_azi + dir.z * sin_azi)/b;
    y = mu*dir.y - a*b*cos_azi;
    z = mu*dir.z + a*(dir.y * dir.z * cos_azi - dir.x * sin_azi)/b;
  }
  return point( x, y, z );
}

