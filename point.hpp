#ifndef _UTIL_POINT_HEADER_
#define _UTIL_POINT_HEADER_

namespace util {

  class point {
    public:
      double x, y, z;

      point( const double x = 0.0, const double y = 0.0, const double z = 0.0 ) :
        x(x), y(y), z(z) {};

      point operator+( const point& p ) const;
      point operator-( const point& p ) const;
      point operator*( const double a ) const;
      point operator/( const double a ) const;
  };

  double dot_product( const point& p, const point& q );
}

#endif
