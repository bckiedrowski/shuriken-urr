#ifndef _MC_MATERIAL_HEADER_
#define _MC_MATERIAL_HEADER_

#include <vector>
#include <map>
#include <memory>

#include "neutronData.hpp"
#include "random.hpp"
#include "xsBuffer.hpp"

namespace mc {

  // material consists of a vector of pairs of atomic densities and xs tables
  class material {
    private:
      uint32_t _size = 0;
      std::vector< uint32_t > _id;
      std::vector< double >   _aden;
      std::vector< std::shared_ptr< xs::neutronData > > _xs;
    public:

      void add( const double aden, std::shared_ptr< xs::neutronData > xs );

      uint32_t size() const { return _size; }
      uint32_t id( const uint32_t k ) const { return _id[k]; }

      double aden( const uint32_t k ) const { return _aden[k]; }
      std::shared_ptr< xs::neutronData > xs( const uint32_t k ) const { return _xs[k]; }

      void fill_total_xs(     const double erg, xs::xsBuffer& xs_buffer, util::RNG::ptr rng ) const;
      void fill_nufission_xs( const double erg, xs::xsBuffer& xs_buffer, util::RNG::ptr rng ) const;

      std::shared_ptr< xs::neutronData > select_isotope(
        const double total_macro, std::vector<double> partial_mirco, util::RNG::ptr rng ) const;
  };

  const auto mat_void = std::make_shared< material >();

}

#endif
