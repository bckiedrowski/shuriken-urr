#include <cmath>

#include "particle.hpp"

void mc::particle::move( const double s ) {
  pos = pos + dir*s;
}

//give particle small nudge along u with a sign given by the dot product
//of the provided vector and the direction of flight
//the intent of this routine is to push a particle into the next cell
//after hitting a surface
void mc::particle::nudge( const util::point u ) {
  const auto s = std::copysign( 1.0e-8, util::dot_product( dir, u ) );
  pos = pos + u*s;
}
