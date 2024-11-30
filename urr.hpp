#ifndef _XS_URR_HEADER_
#define _XS_URR_HEADER_

#include <cassert>
#include <vector>

#include "ace.hpp"
#include "endfMT.hpp"
#include "interpGrid.hpp"

namespace xs {

  struct urrXS {
    uint32_t energy_index, left_xs_index, right_xs_index;
    double total, elastic, fission, ngamma, heating, other;
    double left_elastic_xs,  left_fission_xs,  left_ngamma_xs;
    double right_elastic_xs, right_fission_xs, right_ngamma_xs;
    double interp_factor;
  };

  class urrData {
    private:
      uint32_t _nerg;
      uint32_t _nprob;
      interp_t _interp;
      int32_t  _inelastic_flag;
      int32_t  _other_abs_flag;
      int32_t  _factors_flag;

      //if factors are stored, need to have ownership of a portion of the mean
      //values of the standard cross section table in the unresolved range
      //inelastic and other absorptions reactions may also be stored here if
      //the respective flags are set
      struct background_xs_data {
        std::vector<double> egrid;
        std::vector<double> elastic;
        std::vector<double> fission;
        std::vector<double> ngamma;
        std::vector<double> heating;
        std::vector<double> inelastic;
        std::vector<double> other_abs;
      };
      background_xs_data background_xs;

      std::vector<double> _egrid;
      std::vector< std::vector<double> > _pdf;
      std::vector< std::vector<double> > _cdf;
      std::vector< std::vector<double> > _total;
      std::vector< std::vector<double> > _elastic;
      std::vector< std::vector<double> > _fission;
      std::vector< std::vector<double> > _ngamma;
      std::vector< std::vector<double> > _heating;

      //add specified reaction xs at mt offset to array in xs, which is over
      //the global energy grid segmented over the urr range
      //returns bool indicating whether the reaction was found
      bool tabulate_background( const aceData& acefile, const int32_t start_urr, const int32_t num_erg_urr,
                                const endf_mt_t target_mt, std::vector<double>& xs ) const;
    public:
      struct union_grid {
        std::vector<uint32_t> i1, i2;
        std::vector<double>   pdf, cdf;
      };

      urrData( const aceData& acefile );

      urrXS xs( const double erg, const double random_number ) const;

      double lower_bound() const { return _egrid[0]; }
      double upper_bound() const { return _egrid[_nerg-1]; }

      std::vector<double> egrid() const { return _egrid; };

      uint32_t egrid_size() const { return _nerg; }
      uint32_t num_prob()   const { return _nprob; }

      bool within_urr_range( const double erg ) const {
        return erg >= lower_bound() && erg < upper_bound(); };

      void test();

      std::pair<double,double> table_value( const endf_mt_t rxn, const uint32_t energy_index,
                                            const uint32_t left_index, const uint32_t right_index ) const;

      void adjust( const endf_mt_t mt, const double factor );

      std::vector< union_grid > unionize();
      void analytical_benchmark_1( const double N, const double L );
      void analytical_benchmark_2( const double N, const double L );
      void analytical_benchmark_3( const double N, const double L );
  };

}

#endif
