#ifndef _MC_FIXED_SOURCE_CALC_HEADER_
#define _MC_FIXED_SOURCE_CALC_HEADER_

#include <algorithm>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "cell.hpp"
#include "monteCarloCalc.hpp"
#include "particle.hpp"
#include "point.hpp"

namespace mc {

  class fixedSourceCalc : public monteCarloCalc {
    private:
      std::function< mc::particle( util::RNG::ptr rng ) > _src_dist;

    public:
      fixedSourceCalc( const size_t batch_size, const size_t num_batches,
                       std::function< mc::particle( util::RNG::ptr rng ) > src_dist,
                       const std::vector< std::shared_ptr< mc::estimator > >& estimators,
                       const std::optional< std::shared_ptr< mc::dosEstimator > >   dos_estimators    = std::nullopt,
                       const std::optional< std::shared_ptr< mc::dosUREstimator > > dos_ur_estimators = std::nullopt )
                     : _src_dist(src_dist), monteCarloCalc( batch_size, num_batches, estimators, dos_estimators, dos_ur_estimators ) {};

      bool    is_scoring_batch() const { return true; };
      double  normalization() const {
        return static_cast<double>( _batch_size ) * _batch; };

      particle source_particle( const uint32_t id, util::RNG::ptr rng ) const;

      void score( const estimator_type_t type, const double contribution, const uint32_t geom_id,
                  const mc::particle& p, const xs::xsBuffer& xs_buffer, util::RNG::ptr rng );

      void end_history();
      void end_batch();
  };

}

#endif
