#ifndef _XS_NEUTRON_DATA_HEADER_
#define _XS_NEUTRON_DATA_HEADER_

#include <array>
#include <vector>
#include <limits>
#include <map>
#include <memory>

#include "ace.hpp"
#include "angularDistribution.hpp"
#include "endfMT.hpp"
#include "fissionPhysics.hpp"
#include "point.hpp"
#include "random.hpp"
#include "scatterPhysics.hpp"
#include "urr.hpp"
#include "xsBuffer.hpp"
#include "yield.hpp"

namespace xs {

  class neutronData {
    private:
      struct rxn_t {
        uint32_t ierg, nerg;
        ref_frame_t frame;
        double      qvalue;
        uint32_t    mult;
        std::shared_ptr<yieldData> yield;
        std::shared_ptr<angularDistributionTable> ang_dist;
        std::shared_ptr<scatterPhysics> scatter_physics;
        std::vector<double> xs;
      };

      uint32_t _id;

      uint32_t _zaid;

      double   _temp; //temperature in Kelvin

      uint32_t _nerg;
      uint32_t _nrxn;

      double _emax = 20.0 - 1.0e-11;   //default 20 MeV

      static constexpr int32_t nhash = 4096;   //2^12
      double _hash_min;
      double _hash_max;
      double _hash_delta;
      std::array< int32_t, nhash+1 > _hash_lookup;

      std::vector<double>    _egrid;
      std::vector<endf_mt_t> _mt;

      std::map< endf_mt_t, rxn_t > _rxn;
      std::map< endf_mt_t, rxn_t >::iterator _rxn_total_iter;
      std::map< endf_mt_t, rxn_t >::iterator _rxn_elastic_iter;
      std::map< endf_mt_t, rxn_t >::iterator _rxn_ngamma_iter;
      std::map< endf_mt_t, rxn_t >::iterator _rxn_terminate_iter;

      //note fission given in one of two formats
      //(1) MT = 18, overall fission, store only fission_iter and not the 2nd, 3rd, 4th chances
      //(2) MT = 19,20,21,38 as chances, first chance in fission_iter and others in fission_chances_iter
      std::map< endf_mt_t, rxn_t >::iterator _rxn_fission_iter;
      std::array< std::map< endf_mt_t, rxn_t >::iterator, 3 > _rxn_fission_chance_iter;

      std::vector< std::map< endf_mt_t, rxn_t >::iterator > _other_rxn_iter;

      bool _fissionable = false;
      bool _fission_chances_given = false;

      std::shared_ptr< urrData > _urr;
      bool _urr_exist = false;

      bool read_ace( const aceData& acefile );

      double fast_xs_lookup( const rxn_t& rxn, const uint32_t k, const double r ) const;

      double erg_hash_fxn( const double erg ) const;
    public:
      neutronData( const aceData& acefile, const uint32_t id ) : _id(id) { read_ace( acefile ); };

      uint32_t id()   const { return _id; };
      uint32_t zaid() const { return _zaid; };

      //return Z (charge number), A (mass number), S (state number)
      std::tuple<uint32_t,uint32_t,uint32_t> split_zaid() const;

      void fill_buffer( const double erg, xsBuffer& xs_buffer, util::RNG::ptr rng ) const;

      double find_xs( const endf_mt_t mt, const double erg ) const;
      double yield(   const endf_mt_t mt, const double erg ) const;

      //unresolved resonance functions
      bool urr() const { return  _urr.get(); }

      bool within_urr( const double erg ) const;

      std::shared_ptr< urrData > urr_data() const { return _urr; };

      std::pair<uint32_t,uint32_t> urr_size() const {
        return std::make_pair( _urr->egrid_size(), _urr->num_prob() );
      };
      std::pair<double,double> urr_bounds() const {
        return std::make_pair( _urr->lower_bound(), _urr->upper_bound() ); }
      urrXS  sample_urr( const double erg, const double random_number ) const;
      std::pair<double,double> urr_table( const endf_mt_t rxn, const uint32_t energy_index,
                                          const uint32_t left_index, const uint32_t right_index ) const {
        return _urr->table_value( rxn, energy_index, left_index, right_index ); };

      //random sampling routines
      xs::endf_mt_t select_rxn( const double erg, const xsBuffer& xs_buffer, util::RNG::ptr rng ) const;

      std::pair<double,util::point> sample_secondary(
        const endf_mt_t mt, const double erg, const util::point& dir, util::RNG::ptr rng ) const;

      std::pair<double,util::point> sample_fission(
        const double erg, const util::point& dir, const xsBuffer& xs_buffer, util::RNG::ptr rng ) const;

      //adjustment routine (for testing purposes only)
      void adjust( const endf_mt_t mt, const double factor, const double emin = 0.0,
                   const double emax = std::numeric_limits<double>::max() );
  };

}

#endif
