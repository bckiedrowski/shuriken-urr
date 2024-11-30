#ifndef _MC_ESTIMATOR_HEADER_
#define _MC_ESTIMATOR_HEADER_

#include <vector>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>

#include "endfMT.hpp"
#include "cell.hpp"
#include "xsBuffer.hpp"

namespace mc {

  enum class estimator_type_t {
    SURFACE_CURRENT   = 1,
    TRACK_LENGTH_FLUX = 2,
    COLLISION_FLUX    = 3
  };

  //TODO: energy binning, nu-fission, heating
  class estimator {
    public:
      using ptr = std::shared_ptr< estimator >;
    private:
      std::string      _name;
      int32_t          _id;
      estimator_type_t _type;

      std::set<uint32_t> _geom_id;                         //=empty means all cells/surfaces
      uint32_t           _zaid    = 0;                     //=zero denotes all isotopes
      xs::endf_mt_t      _mt      = xs::endf_mt_t::NONE;   //=default of no multiplier for flux/current

      uint32_t _nebins = 0;
      std::vector< double > _egrid;

      bool _within_xs_buffer = false;

      std::vector<double> _sum_hist;   //running total for the history

      std::vector<double> _mean;   //sum of history scores
      std::vector<double> _stdv;   //sum square of history scores
    public:
      estimator( const std::string name, const int32_t id, const estimator_type_t type,
                 const std::optional< std::set<uint32_t> > geom_id,
                 const std::optional< uint32_t >           zaid,
                 const std::optional< xs::endf_mt_t >      mt,
                 const std::optional< std::vector<double> >& ebins );

      std::string      name() const { return _name; }
      int32_t          id()   const { return _id;   }
      estimator_type_t type() const { return _type; }

      uint32_t      zaid()    const { return _zaid; }
      xs::endf_mt_t mt()      const { return _mt; }
      uint32_t nebins()       const { return _nebins; }

      bool within_geom_set( const uint32_t geom_id ) const;
      bool within_energy_range( const double erg )   const;

      //returns energy bin index (1 = lowest energy bin such that total is indexed by 0)
      uint32_t energy_bin( const double erg ) const;

      double egrid_value( const uint32_t k ) const { return _egrid[k]; };

      void score( const double contribution, const uint32_t geom_id, const double erg,
                  const std::shared_ptr< geom::cell > cell_ptr,
                  const xs::xsBuffer& xs_buffer );

      void end_history();

      //returns current accumulated score for this history for the specified bin
      double current_score( const uint32_t k ) const { return _sum_hist[k]; }

      double mean( const double normalization, const uint32_t k ) const;
      double stdv( const double normalization, const uint32_t k ) const;
      void   write( const double normalization ) const;
  };

  void score_estimators_of_type( const estimator_type_t type,
    const std::vector< std::shared_ptr<estimator> >& estimators,
    const double contribution, const uint32_t geom_id, const double erg,
    const std::shared_ptr< geom::cell > cell_ptr, const xs::xsBuffer& xs_buffer );

  void end_history_estimators( const std::vector< std::shared_ptr<estimator> >& estimators );

  void write_estimators( const std::vector< std::shared_ptr<estimator> >& estimators, const double normalization );

}
#endif
