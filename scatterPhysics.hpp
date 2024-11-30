#ifndef _XS_SCATTER_PHYSICS_HEADER_
#define _XS_SCATTER_PHYSICS_HEADER_

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

namespace xs {

  enum class ref_frame_t { UNKNOWN, COM_FRAME, LAB_FRAME };

  class scatterPhysics {
    protected:
      bool _present = false;
      bool _elastic = false;

      ref_frame_t _ref_frame;

      //reduced version of the interpGrid class used in ACE laws 1,4,44,61 to
      //interpolate between incident energy grids
      //supports histogram and lin-lin interpolation with method to return
      //lower energy grid index and the interpolation fraction
      class energyGridInterp {
        private:
          std::vector<uint32_t> interp_bounds;
          std::vector<interp_t> interp_types;
          std::vector<double>   egrid;         //grid incident energies
        public:
          energyGridInterp( const aceData& acefile, const uint32_t loc );

          //returns incident energy grid (lower) index and interpolation fraction
          std::pair< uint32_t, double > search( const double erg ) const;
      };

      class aceLaw {
        protected:
          ref_frame_t ref_frame;
          double A;   //= atomic weight ratio
        public:
          virtual std::pair<double,util::point> sample(
            const double erg, const util::point& dir, util::RNG::ptr rng ) = 0;

          std::pair<double,double> com_to_lab( const double erg_in, const double erg_cm,
                                               const double mu_cm,  const double awr     ) const;

      };

      class freeGasLaw : public aceLaw {
        private:
          double temp; //table temperature

          std::shared_ptr<angularDistributionTable> ang_dist;
        public:
          freeGasLaw( const aceData& acefile, std::shared_ptr<angularDistributionTable> ang );

          std::pair<double,util::point> sample( const double erg, const util::point& dir, util::RNG::ptr rng );
      };

      class aceLaw3 : public aceLaw {
        private:
          double t1;  //= (A+1)/A * abs(Q) = LDIS(1)
          double t2;  //= (A/(A+1))^2      = LDIS(2)

          std::shared_ptr<angularDistributionTable> ang_dist;
        public:
          aceLaw3( const aceData& acefile, const uint32_t loc, ref_frame_t frame,
                   std::shared_ptr<angularDistributionTable> ang );

          std::pair<double,util::point> sample( const double erg, const util::point& dir, util::RNG::ptr rng );
      };

      //class for ACE laws 4,44,61
      // 4 = tabular for energy-angle sampled independently
      //44 = tabular for angle sampled with kallbach-mann systematics correlated with energy
      //61 = tabular for angle sampled with arbitrary tabular distributions correlated with energy
      //note: could be done with inheritence, but would be more complicated given the different
      //      ways the data are represented and angles are sampled
      enum class ace_tabular_t { LAW4, LAW44, LAW61 };
      class aceLawTabular : public aceLaw {
        private:
          ace_tabular_t law_type;
          uint32_t ndiscrete = 0;

          struct outgoing_data {
            interp_t interp;
            std::vector<double> eout, pdf, cdf;
            //angular data for kallbach-mann (law44)
            std::vector<double> A, R;
            //angular data for tabular distributions (law61)
            std::vector< std::shared_ptr<angDist> > ang_dist;
          };
          std::vector<uint32_t> locators;
          std::vector<outgoing_data> edata; //outgoing energy data

          energyGridInterp inc_egrid;

          //assigned for law4, unassigned (nullptr) for law44 and law61
          std::shared_ptr<angularDistributionTable> ang_dist;
        public:
          aceLawTabular( ace_tabular_t type, const aceData& acefile, const uint32_t loc_dlw,
            const uint32_t loc_ldat, ref_frame_t frame, std::shared_ptr<angularDistributionTable> ang = nullptr );

          std::pair<double,util::point> sample( const double erg, const util::point& dir, util::RNG::ptr rng );
      };

      class aceLaw7 : public aceLaw {
        private:
          double      restrict_energy;
          interpGrid  interp;

          std::shared_ptr<angularDistributionTable> ang_dist;
        public:
          aceLaw7( const aceData& acefile, const uint32_t loc, ref_frame_t frame,
                   std::shared_ptr<angularDistributionTable> ang );

          std::pair<double,util::point> sample( const double erg, const util::point& dir, util::RNG::ptr rng );
      };

      class aceLaw9 : public aceLaw {
        private:
          double      restrict_energy;
          interpGrid  interp;

          std::shared_ptr<angularDistributionTable> ang_dist;
        public:
          aceLaw9( const aceData& acefile, const uint32_t loc, ref_frame_t frame,
                   std::shared_ptr<angularDistributionTable> ang );

          std::pair<double,util::point> sample( const double erg, const util::point& dir, util::RNG::ptr rng );
      };

      struct lawBlock {
        uint32_t   id = 0;
        std::shared_ptr<interpGrid> prob;
        std::shared_ptr<aceLaw>     data;
      };

      std::shared_ptr<angularDistributionTable> _ang_dist;
      std::vector<lawBlock> _laws;

      void create_laws( const aceData& acefile, const uint32_t loc_dlwp,
                        const uint32_t loc_dlw, const uint32_t mt_offset );
    public:
      scatterPhysics( const aceData& acefile, const uint32_t mt_offset );

      bool operator()() const { return _present; }

      //return energy and direction vector (lab frame) pair
      virtual std::pair<double,util::point> sample( const double erg, const util::point& dir, util::RNG::ptr rng ) const;
  };
}

#endif
