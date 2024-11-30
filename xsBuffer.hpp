#ifndef _XS_XS_BUFFER_HEADER_
#define _XS_XS_BUFFER_HEADER_

#include <vector>
#include <memory>

#include "endfMT.hpp"
#include "random.hpp"

namespace xs {

  class neutronData;

  struct urr_cache_t {
    uint32_t energy_index, left_xs_index, right_xs_index;
    double   interp_factor, left_elastic_xs,  left_fission_xs,  left_ngamma_xs,  left_nufission_xs,
                            right_elastic_xs, right_fission_xs, right_ngamma_xs, right_nufission_xs;
  };

  //note vectors ordered by global unique isotope ID number
  class xsBuffer {
    public:
      xsBuffer( const uint32_t niso );

      //arrays for micoscopic cross sections for each isotope in global problem
      std::vector<double> erg, total_xs, nufission_xs, elastic_xs, fission_xs, ngamma_xs, urr_rn;

      //macroscopic cross section data specific to the current cell
      double total_macro_xs     = 0.0;
      double nufission_macro_xs = 0.0;

      //isotope zaid and reaction type for collision
      uint32_t  collided_zaid = 0;
      endf_mt_t collided_mt   = endf_mt_t::NONE;

      double operator()( const xs::endf_mt_t mt, const uint32_t id ) const;

      std::vector<urr_cache_t> urr_data;
      uint32_t urr_energy_index(   const uint32_t id ) const { return urr_data[id].energy_index;   };
      uint32_t urr_left_xs_index(  const uint32_t id ) const { return urr_data[id].left_xs_index;  };
      uint32_t urr_right_xs_index( const uint32_t id ) const { return urr_data[id].right_xs_index; };
      double   urr_interp_factor(  const uint32_t id ) const { return urr_data[id].interp_factor;  };

      std::pair< double, double > urr_table_xs( const xs::endf_mt_t mt, const uint32_t id ) const;

      void flush();
  };

}

#endif
