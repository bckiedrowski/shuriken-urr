#ifndef _XS_FISSION_PHYSICS_HEADER_
#define _XS_FISSION_PHYSICS_HEADER_

#include <array>
#include <vector>
#include <memory>

#include "ace.hpp"
#include "angularDistribution.hpp"
#include "constants.hpp"
#include "endfMT.hpp"
#include "interpGrid.hpp"
#include "point.hpp"
#include "random.hpp"
#include "scatterPhysics.hpp"
#include "yield.hpp"

namespace xs {

  class fissionPhysics : public scatterPhysics {
    private:
      //number of delayed precursor groups
      uint32_t _num_delayed_groups = 0;

      //ratio of prompt to total yields (nu-bar) give the probability of prompt fission
      std::shared_ptr< yieldData > _total_yield, _prompt_yield;

      //data in _delayed_group_prob are probabilities denoting the probability that a neutron
      //is emitted from a precusor group conditional on the fission being delayed
      //this implies the sum of the probabilities should be one
      std::vector< std::shared_ptr<interpGrid> > _delayed_group_prob;

      //decay constants are currently unused but could easily advance the particle state in time
      //but would require changing the return type of the sample routines (probably particle)
      //note ace file gives these in per shake, here convert to per second
      std::vector< double > _decay_constants;

      //spectra for both prompt and delayed neutrons are sampled from inelastic scattering laws
      //law(s) for prompt fission are in the _laws vector of the base scatterPhysics class first
      //the _group_indices vector gives offsets into the law array of valid laws to search
      //for delayed spectra (each of which could have their own set of laws)
      //size is number of groups + 2 and structured as follows:
      //[0] = location of prompt fission laws (normally zero)
      //[1] through [_num_delayed_groups] = location of delayed fission laws
      //[_num_delayed_groups+1] = size of the array
      //this ordering makes the loop for the law selection simpler
      std::vector< uint32_t > _group_law_offsets;

    public:
      fissionPhysics( const aceData& acefile, const uint32_t mt_offset );

      std::pair<double,util::point> sample( const double erg, const util::point& dir, util::RNG::ptr rng ) const override;
  };
}

#endif
