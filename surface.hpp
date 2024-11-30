#ifndef _GEOM_SURFACE_HEADER_
#define _GEOM_SURFACE_HEADER_

#include <cstdint>
#include <string>

#include "point.hpp"

namespace geom {

  enum class surface_bc_t {
    NONE       = 0,
    REFLECTING = 1
  };

  class surface {
    private:
      std::string  _name;
      uint32_t     _id;
      surface_bc_t _bc;

      virtual util::point grad( const util::point p ) const = 0;
    public:
      surface( const std::string name, const uint32_t id, const surface_bc_t bc = surface_bc_t::NONE ) :
        _name(name), _id(id), _bc(bc) {};

      std::string  name() const { return _name; }
      uint32_t     id()   const { return _id;   }
      surface_bc_t bc()   const { return _bc;   }

      virtual double eval( const util::point p ) const = 0;
      virtual double distance( const util::point p, const util::point u ) const = 0;

      util::point normal(  const util::point p ) const;
      util::point reflect( const util::point p, const util::point u ) const;
  };

  class plane : public surface {
    private:
      double a, b, c, d;
      util::point grad( const util::point p ) const override;
    public:
      plane( const std::string name, const uint32_t id, const double a, const double b,
             const double c, const double d, const surface_bc_t bc = surface_bc_t::NONE )
           : surface(name,id,bc), a(a), b(b), c(c), d(d) {};

      double eval( const util::point p ) const override;
      double distance( const util::point p, const util::point u ) const override;
  };

  class sphere : public surface {
    private:
      double x0, y0, z0, r;
      util::point grad( const util::point p ) const override;
    public:
      sphere( const std::string name, const uint32_t id, const double x0, const double y0,
              const double z0, const double r, const surface_bc_t bc = surface_bc_t::NONE )
            : surface(name,id,bc), x0(x0), y0(y0), z0(z0), r(r) {};

      double eval( const util::point p ) const override;
      double distance( const util::point p, const util::point u ) const override;
  };

  class cylinder_z : public surface {
    private:
      double x0, y0, r;
      util::point grad( const util::point p ) const override;
    public:
      cylinder_z( const std::string name, const uint32_t id, const double x0, const double y0,
                  const double r, const surface_bc_t bc = surface_bc_t::NONE )
                : surface(name,id,bc), x0(x0), y0(y0), r(r) {};

      double eval( const util::point p ) const override;
      double distance( const util::point p, const util::point u ) const override;
  };

}

#endif
