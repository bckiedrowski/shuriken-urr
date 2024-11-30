#ifndef _GEOM_CELL_HEADER_
#define _GEOM_CELL_HEADER_

#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "material.hpp"
#include "point.hpp"
#include "surface.hpp"

namespace geom {

  struct nextSurface {
    double distance = std::numeric_limits<double>::max();
    std::shared_ptr< surface > surface = nullptr;
  };

  //assume pure intersection operator for surfaces for now
  class cell {
    private:
      std::string _name;
      uint32_t    _id;
      std::vector< std::shared_ptr< surface > > _surfaces;
      std::vector< int > _senses;

      std::shared_ptr< mc::material > _material;
    public:
      cell( const std::string name, const uint32_t id, const std::vector< std::shared_ptr< surface > > surfaces,
            const std::vector< int > senses, std::shared_ptr<mc::material> material ) :
            _name(name), _id(id), _surfaces(surfaces), _senses(senses), _material(material) {};

      std::string name() const { return _name; }
      uint32_t    id()   const { return _id;   }
      std::shared_ptr< mc::material > material() const { return _material; }

      bool within( util::point p ) const;
      nextSurface next_surface( util::point p, util::point u ) const;
  };

  std::shared_ptr<cell> find_cell( const util::point p, const std::vector<std::shared_ptr<cell>> cells );

}

#endif
