#ifndef _MC_MONTE_CARLO_CALC_HEADER_
#define _MC_MONTE_CARLO_CALC_HEADER_

#include <algorithm>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "calcMode.hpp"
#include "cell.hpp"
#include "dosEstimator.hpp"
#include "dosUREstimator.hpp"
#include "estimator.hpp"
#include "particle.hpp"
#include "point.hpp"

namespace mc {

  class monteCarloCalc {
    protected:
      size_t _batch = 0;
      size_t _batch_size;
      size_t _num_batches;

      std::vector< std::shared_ptr< mc::estimator > > _estimators;
      std::shared_ptr< mc::dosEstimator   > _dos_estimators;
      std::shared_ptr< mc::dosUREstimator > _dos_ur_estimators;
    public:
      monteCarloCalc( const size_t batch_size, const size_t num_batches,
                      const std::vector< std::shared_ptr< mc::estimator > >& estimators,
                      const std::optional< std::shared_ptr< mc::dosEstimator > >   dos_estimators    = std::nullopt,
                      const std::optional< std::shared_ptr< mc::dosUREstimator > > dos_ur_estimators = std::nullopt )
                    : _batch_size(batch_size), _num_batches(num_batches), _estimators(estimators),
                      _dos_estimators( dos_estimators.value_or( nullptr ) ),
                      _dos_ur_estimators( dos_ur_estimators.value_or( nullptr ) ) {};

      size_t  batch_size()  const { return _batch_size;  };
      size_t  num_batches() const { return _num_batches; };
      bool    is_done()     const { return _batch >= _num_batches; };

      virtual bool   is_scoring_batch() const = 0;
      virtual double normalization()    const = 0;

      virtual particle source_particle( const uint32_t id, util::RNG::ptr rng ) const = 0;

      virtual void score( const estimator_type_t type, const double contribution, const uint32_t geom_id,
                          const mc::particle& p, const xs::xsBuffer& xs_buffer, util::RNG::ptr rng ) = 0;

      virtual void end_history() = 0;
      virtual void end_batch()   = 0;
  };

}

#endif
