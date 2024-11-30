#ifndef _XS_ANGULAR_DISTRIBUTION_HEADER_
#define _XS_ANGULAR_DISTRIBUTION_HEADER_

#include <array>
#include <vector>
#include <memory>

#include "ace.hpp"
#include "endfMT.hpp"
#include "interpGrid.hpp"
#include "random.hpp"

namespace xs {

  class angDist {
    public:
      virtual double sample( util::RNG::ptr rng ) const = 0;
  };

  class isotropicAngDist : public angDist {
    public:
      isotropicAngDist() {};
      double sample( util::RNG::ptr rng ) const;
  };

  class equiprobable32BinAngDist : public angDist {
    private:
      std::array<double, 33> _mu;
    public:
      equiprobable32BinAngDist( const aceData& acefile, const uint32_t loc );
      double sample( util::RNG::ptr rng ) const;
  };

  class tabularAngDist : public angDist {
    private:
      interp_t _interp_type;
      uint32_t _ncos;
      std::vector<double> _mu;
      std::vector<double> _pdf;
      std::vector<double> _cdf;
    public:
      tabularAngDist( const aceData& acefile, const uint32_t loc );
      double sample( util::RNG::ptr rng ) const;
  };

  class angularDistributionTable {
    private:
      bool     _present   = false;
      bool     _isotropic = false;
      uint32_t _nerg;
      std::vector<double> _egrid;
      std::vector< std::shared_ptr<angDist> > _ang;

    public:
      angularDistributionTable( const aceData& acefile, const uint32_t mt_offset );

      bool operator()() const { return _present; }
      double sample( const double erg, util::RNG::ptr rng ) const;
  };
}

#endif
