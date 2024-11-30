#ifndef _MC_PARTICLE_HEADER_
#define _MC_PARTICLE_HEADER_

#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "point.hpp"
#include "cell.hpp"

namespace mc {

  class particle {
    private:
      bool _alive = true;
    public:
      util::point pos, dir;
      double erg;
      double wgt = 1.0;

      std::shared_ptr< geom::cell > cell;

      bool alive() const { return _alive; }

      void move( const double s );
      void nudge( const util::point u );

      void kill() { _alive = false; };
  };

}

#endif
