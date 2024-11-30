#ifndef _MC_K_EIGENVALUE_CALC_HEADER_
#define _MC_K_EIGENVALUE_CALC_HEADER_

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "cell.hpp"
#include "monteCarloCalc.hpp"
#include "particle.hpp"
#include "point.hpp"

namespace mc {

  struct fissionSite {
    util::point pos;
    double      erg;
    std::shared_ptr< geom::cell > cell_ptr;

    uint32_t    dos_id = 0;

    fissionSite() {};
    fissionSite( const util::point pos, const double erg,
                 const std::shared_ptr< geom::cell > cell_ptr ) :
      pos(pos), erg(erg), cell_ptr(cell_ptr) {};
  };

  class kEigenvalueCalc : public monteCarloCalc {
    private:
      double _keff_norm = 1.0;
      double _mean_keff = 0.0;
      double _stdv_keff = 0.0;
      double _keff_est  = 0.0;
      double _source_weight = 1.0;

      size_t _batch_size_initial;
      size_t _num_inactive;

      size_t _nbanked = 0;

      std::vector< mc::fissionSite > _source_bank;
      std::vector< mc::fissionSite > _fission_bank;
    public:
      kEigenvalueCalc( const size_t batch_size, const size_t num_inactive, const size_t num_batches,
                       const util::point pos, const std::vector< std::shared_ptr<geom::cell> >& cells,
                       const std::vector< std::shared_ptr< mc::estimator > >& estimators,
                       const std::optional< std::shared_ptr< mc::dosEstimator > >   dos_estimators    = std::nullopt,
                       const std::optional< std::shared_ptr< mc::dosUREstimator > > dos_ur_estimators = std::nullopt );

      bool    is_scoring_batch() const { return _batch >= _num_inactive; };
      double  normalization() const {
        return static_cast<double>( _batch_size_initial ) * std::max( 0, static_cast<int>(_batch) - static_cast<int>(_num_inactive) ); };

      particle source_particle( const uint32_t id, util::RNG::ptr rng ) const;

      void score( const estimator_type_t type, const double contribution, const uint32_t geom_id,
                  const mc::particle& p, const xs::xsBuffer& xs_buffer, util::RNG::ptr rng );

      void end_history();
      void end_batch();
  };

}

#endif
