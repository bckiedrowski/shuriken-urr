#ifndef _MC_DOS_UR_ESTIMATOR_HEADER_
#define _MC_DOS_UR_ESTIMATOR_HEADER_

#include <vector>
#include <map>
#include <memory>
#include <optional>

#include "dosEstimator.hpp"
#include "endfMT.hpp"
#include "estimator.hpp"
#include "cell.hpp"
#include "xsBuffer.hpp"

namespace mc {

  //differential operator sampling perturbation estimator class for sensitivity coefficients
  class dosUREstimator {
    public:
      using rxn_ratio_vector = std::vector< std::pair< mc::estimator::ptr, mc::estimator::ptr > >;
    private:
      //calculation mode
      mc::mc_calc_t _calc_mode;

      //peturbation parameters
      std::vector< std::string   >        _name;
      std::vector< std::set< uint32_t > > _cell_id;
      std::vector< uint32_t >             _zaid;
      std::vector< xs::endf_mt_t >        _mt;

      std::vector< uint32_t > _nebins;
      std::vector< uint32_t > _nprob;
      std::vector< std::vector< double > > _egrid;

      uint32_t _nprofiles   = 0;       //=number of sensitivity profiles per response
      uint32_t _nparameters = 0;       //=total number of sensitivity coefficients per response

      std::vector< uint32_t > _parameter_id;  //=map from profile index to parameter index

      //response data
      static constexpr int32_t _internal_keff_estimator_id = -1;  //keigenvalue mode has internal keff estimator
      uint32_t _nestimators = 0;       //=number of unique user-defined estimators for sensitivities
      uint32_t _nresponses  = 0;       //=number of unique response sensitivities
      uint32_t _nrxn_ratios = 0;       //=number of user-defined reaction ratios

      std::map< int32_t, mc::estimator::ptr >  _estimators;
      rxn_ratio_vector _rxn_ratios;

      std::vector< uint32_t > _response_id;  //=map to response index in order that estimators were placed

      std::map< int32_t, std::vector< std::vector< double > > > _sens_mean;
      std::map< int32_t, std::vector< std::vector< double > > > _sens_stdv;

      //TODO: integrate mean with _responses
      //double _keff_sens_mean = 0.0;
      //double _keff_sens_stdv = 0.0;

      //cache perturbed reaction cross section
      //computed upon stream and may be needed again at collision
      std::vector< double > _pert_rxn_xs_cache;       //TODO: probably delete...
      std::vector< double > _pert_rxn_left_xs_cache;
      std::vector< double > _pert_rxn_right_xs_cache;

      //cache energy index offset to avoid searching multiple times
      std::vector< uint32_t > _erg_offset_cache;
      std::vector< uint32_t > _iso_id_cache;

      //flag for whether to check the cross section buffer as in the case of total, fission,
      //elastic, or (n,gamma) cross sections should they be perturbed, otherwise need to do
      //a more expensive lookup on the cross section table
      //note that this is essentially required to handle unresolved resonances correctly
      std::vector< bool > _within_xs_buffer;

      //indirect and response direct effects for current history
      std::vector< double > _indirect_effect;
      std::vector< std::vector< double > > _direct_effect;

      //fission source renormalization
      std::vector< double > _total_perturbed_weight;
      std::vector< double > _fission_source_renorm;

      //vector of derivative estimates stored for each fission site
      std::vector< std::vector< double > > _fission_bank;
      std::vector< std::vector< double > > _source_bank;

      //counter for fission bank index returned from insert_fission_source for storage
      //within the source bank
      uint32_t _fbank_index = 0;

    public:
      dosUREstimator( const std::map< uint32_t, std::shared_ptr<xs::neutronData> >& xs_map,
                      const std::vector< dos_profile_t > parameters,
                      const std::optional< rxn_ratio_vector > rxn_ratios,
                      const mc::mc_calc_t calc_mode );

      void resize_source_bank( const uint32_t bank_size );

      //ipf = sensitivity profile index
      bool within_energy_range( const double erg, const uint32_t ipf ) const;

      //called to setup the history by pulling the source perturbation weight and
      //normalizing
      void start_history( const uint32_t bank_index );

      //scoring functions
      void score_indirect_stream( const double wgt_times_distance, const double erg,
                                  const std::shared_ptr<geom::cell> cell_ptr,
                                  const xs::xsBuffer& xs_buffer );

      void score_indirect_collision( const double wgt, const double erg,
                                     const std::shared_ptr<geom::cell> cell_ptr,
                                     const xs::xsBuffer& xs_buffer );

      void score_direct_response( const double wgt_times_distance, const double erg,
                                  const std::shared_ptr<geom::cell> cell_ptr,
                                  const xs::xsBuffer& xs_buffer );

      //inserts the derivative estimate, indirect + direct at fission production that is computed
      //by this routine
      uint32_t insert_fission_source( const double wgt, const double erg,
                                      const std::shared_ptr<geom::cell> cell_ptr,
                                      const uint32_t num_neutrons,
                                      const xs::xsBuffer& xs_buffer );

      void end_history();

      //sets structure up for the next cycle: computes normalization, swaps fission and source banks
      //and zeros out tallies
      void end_batch( const uint32_t total_neutrons_banked );

      void write_sensitivity( const double normalization ) const;
  };

}

#endif
