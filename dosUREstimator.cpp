#include <cmath>
#include <iomanip>
#include <iostream>

#include "binarySearch.hpp"

#include "dosUREstimator.hpp"

mc::dosUREstimator::dosUREstimator( const std::map< uint32_t, std::shared_ptr<xs::neutronData> >& xs_map,
                                    const std::vector< dos_profile_t > profiles,
                                    const std::optional< rxn_ratio_vector > rxn_ratios,
                                    const mc::mc_calc_t calc_mode ) {
  _calc_mode = calc_mode;
  _nprofiles   = profiles.size();

  _name.resize( _nprofiles );
  _cell_id.resize( _nprofiles );
  _zaid.resize( _nprofiles );
  _mt.resize( _nprofiles );
  _nebins.resize( _nprofiles );
  _nprob.resize( _nprofiles );
  _egrid.resize( _nprofiles );
  _parameter_id.resize( _nprofiles );
//  _within_xs_buffer.resize( _nprofiles );

  _pert_rxn_xs_cache.resize(       _nprofiles, 0.0 );
  _pert_rxn_left_xs_cache.resize(  _nprofiles, 0.0 );
  _pert_rxn_right_xs_cache.resize( _nprofiles, 0.0 );

  _erg_offset_cache.resize( _nprofiles, 0 );
  _iso_id_cache.resize( _nprofiles, 0 );

  //number of parameters is the sum of the size of each of the profiles
  //a profile is dimensioned by the number of energy grid points times probability points
  _nparameters = 0;
  for ( auto pr : profiles ) {
    if ( ! pr.zaid ) {
      std::cout << " must specify isotope for unresolved resonance sensitivity estimator " << '\n';
      assert( false );
    }
    if ( ! pr.mt ) {
      std::cout << " must specify reaction type for unresolved resonance sensitivity estimator " << '\n';
      assert( false );
    }
    if ( pr.egrid ) {
      std::cout << " energy grid not permitted for unresolved resonance sensitivity estimator " << '\n';
    }
    //find isotope in the cross section map
    const auto zaid  = pr.zaid.value();
    const auto it    = xs_map.find( zaid );
    if ( it == xs_map.end() ) {
      std::cout << " isotope " << std::to_string(zaid) << " not found in problem " << '\n';
      assert( false );
    }
    const auto urdat  = it->second->urr_data();
    const auto nebins = urdat->egrid_size();
    const auto nprob  = urdat->num_prob();

    //each profile has a parameter for:
    //(1) the sum total over the entire table                      = 1
    //(2) the sum over all probabilities at each energy grid point = nebins
    //(3) each cross section or factor on the probability table    = nebins * nprob
    //ordering: { sum total, total over ebin 0, ebin 0/prob 0, ebin 0/ prob 1, ...
    //                       total over ebin 1, ebin 1/prob 0, ebin 1/ prob 2, ... }
    _nparameters += 1 + nebins + nebins * nprob;
  }

  auto ip = 0;   //ip = parameter index, ipf = profile index
  for ( auto ipf = 0 ; ipf < _nprofiles ; ++ipf ) {
    _name[ipf]    = profiles[ipf].name;
    _cell_id[ipf] = profiles[ipf].cell_id.value_or( std::set<uint32_t>() );
    _zaid[ipf]    = profiles[ipf].zaid.value_or( 0 );
    _mt[ipf]      = profiles[ipf].mt.value_or( xs::endf_mt_t::TOTAL );

/*
    //only elastic, fission, and n,gamma acceptable
    if ( _mt[ipf] != xs::endf_mt_t::ELASTIC && _mt[ipf] != xs::endf_mt_t::FISSION
      && _mt[ipf] != xs::endf_mt_t::Z_GAMMA ) {
      std::cout << " only elastic, fission, or (n,gamma) allowed for unresolved resonance sensitivity estimator " << '\n';
      assert( false );
    }
*/
    //only total. elastic, fission, and n,gamma acceptable
    if ( _mt[ipf] != xs::endf_mt_t::ELASTIC && _mt[ipf] != xs::endf_mt_t::FISSION
      && _mt[ipf] != xs::endf_mt_t::Z_GAMMA && _mt[ipf] != xs::endf_mt_t::TOTAL  ) {
      std::cout << " only total, elastic, fission, or (n,gamma) allowed for unresolved resonance sensitivity estimator " << '\n';
      assert( false );
    }

    const auto it    = xs_map.find( _zaid[ipf] );
    const auto urdat = it->second->urr_data();
    _nebins[ipf] = urdat->egrid_size();
    _nprob[ipf]  = urdat->num_prob();

    _parameter_id[ipf] = ip;
    //ip += 1 + ( _nebins[ipf] + 1 ) * _nprob[ipf];
    ip += 1 + _nebins[ipf] + _nebins[ipf] * _nprob[ipf];

    _egrid[ipf] = urdat->egrid();
  }

  //create a special keff response by default (index by -1 in map) for k-eigenvalue calculations
  if ( _calc_mode == mc::mc_calc_t::K_EIGENVALUE ) {
    _estimators[_internal_keff_estimator_id] = std::make_shared< mc::estimator >(
      "keff", _internal_keff_estimator_id, mc::estimator_type_t::TRACK_LENGTH_FLUX,
      std::nullopt, std::nullopt, xs::endf_mt_t::NUFISSION, std::nullopt );

    _sens_mean[_internal_keff_estimator_id] =
      std::vector< std::vector<double> >( 1, std::vector<double>( _nparameters, 0.0 ) );
    _sens_stdv[_internal_keff_estimator_id] =
      std::vector< std::vector<double> >( 1, std::vector<double>( _nparameters, 0.0 ) );

    _nresponses = 1;
  }

  //create user defined responses from ratios vector
  if ( rxn_ratios ) {
    _rxn_ratios = rxn_ratios.value();
    for ( auto rrit : _rxn_ratios ) {
      // check that taking ratio makes sense and meets one of the following
      // (1) numerator and denominator have no energy bins (first condition)
      // (2) numerator has energy bins, but the denominator does not (second condition)
      // (3) numerator and denominator have same number of energy bins (first condition)
      //note that denominator does not need to be specified for fixed-source calculations to make sense
      if (   rrit.second ) assert( ( rrit.first->nebins() == rrit.second->nebins() ) || rrit.second->nebins() == 0 );
      if ( ! rrit.second ) assert( _calc_mode == mc::mc_calc_t::FIXED_SOURCE );

      if ( _estimators.find( rrit.first->id() )  == _estimators.end() ) {
        _estimators[ rrit.first->id() ] = rrit.first;

        const auto nr = rrit.first->nebins() >= 2 ? rrit.first->nebins() + 1 : 1;
        _nresponses += nr;

        _sens_mean[ rrit.first->id() ] =
          std::vector< std::vector<double> >( nr, std::vector<double>( _nparameters, 0.0 ) );
        _sens_stdv[ rrit.first->id() ] =
          std::vector< std::vector<double> >( nr, std::vector<double>( _nparameters, 0.0 ) );
      }
      if ( rrit.second ) {
        //denominator has been specified, if nullptr then we treat as a constant (only valid for fixed source)
        if ( _estimators.find( rrit.second->id() ) == _estimators.end() ) {
          _estimators[ rrit.second->id() ] = rrit.second;

          const auto nr = rrit.second->nebins() >= 2 ? rrit.second->nebins() + 1 : 1;
          _nresponses += nr;

          _sens_mean[ rrit.second->id() ] =
            std::vector< std::vector<double> >( nr, std::vector<double>( _nparameters, 0.0 ) );
          _sens_stdv[ rrit.second->id() ] =
            std::vector< std::vector<double> >( nr, std::vector<double>( _nparameters, 0.0 ) );
        }
      }
    }
  }
  _nestimators = _estimators.size();

  //error if no estimators
  assert( _nestimators > 0 );

  //loop again to ensure consistent ordering
  _response_id.resize( _nestimators, 0 );
  auto ir = 0, iest = 0;
  for ( auto est : _estimators ) {
    _response_id[iest] = ir;
    ir += est.second->nebins() >= 2 ? est.second->nebins() + 1 : 1;
    ++iest;
  }

  _indirect_effect.resize( _nparameters, 0.0 );
  _direct_effect.resize( _nresponses, std::vector<double>( _nparameters, 0.0 ) );
}

bool mc::dosUREstimator::within_energy_range( const double erg, const uint32_t ipf ) const {
  if ( erg < _egrid[ipf][0] || erg >= _egrid[ipf][_nebins[ipf]-1] ) return false;
  return true;
}

void mc::dosUREstimator::resize_source_bank( const uint32_t bank_size ) {
  assert( _calc_mode == mc::mc_calc_t::K_EIGENVALUE );

  _fission_source_renorm.resize( _nparameters, 0.0 );
  _total_perturbed_weight.resize( _nparameters, 0.0 );

  _fission_bank.resize( bank_size, std::vector< double >( _nparameters, 0.0 ) );
  _source_bank.resize(  bank_size, std::vector< double >( _nparameters, 0.0 ) );
}

void mc::dosUREstimator::start_history( const uint32_t bank_index ) {
  //put fission source perturbation (with renormalization) from history that produced
  //fission neutron into the indirect effect vector of parameters
  for ( auto ip = 0 ; ip < _nparameters ; ++ip ) {
    if ( _calc_mode == mc::mc_calc_t::K_EIGENVALUE ) {
      //put fission source perturbation (with renormalization) from history that produced
      //fission neutron into the indirect effect vector of parameters
      _indirect_effect[ip] = _source_bank[bank_index][ip] - _fission_source_renorm[ip];
    }
    else {
      //fixed source assumes no perturbation to source term
      _indirect_effect[ip] = 0.0;
    }
  }
  //zero out direct effect terms for each reponse parameter pair
  //note this can be expensive for a large number of parameters, so important to loop
  //over in the correct order
  for ( auto ir = 0 ; ir < _direct_effect.size() ; ++ir ) {
    for ( auto ip = 0 ; ip < _nparameters ; ++ip ) {
      _direct_effect[ir][ip] = 0.0;
    }
  }
}

void mc::dosUREstimator::score_indirect_stream( const double wgt_times_distance, const double erg,
                                              const std::shared_ptr<geom::cell> cell_ptr,
                                              const xs::xsBuffer& xs_buffer ) {
  //grab material and filter out voids
  //TODO: verify we do not need to zero out cache for void
  const auto mat = cell_ptr->material();
  if ( ! mat ) return;

  for ( auto ipf = 0 ; ipf < _nprofiles ; ++ipf ) {
    const auto ip = _parameter_id[ipf];
    _pert_rxn_xs_cache[ipf] = 0.0;
    _erg_offset_cache[ipf]  = 0;

    _pert_rxn_left_xs_cache[ipf]  = 0.0;
    _pert_rxn_right_xs_cache[ipf] = 0.0;

    //filter out not being in impacted cell
    if ( ! _cell_id[ipf].empty() && _cell_id[ipf].find( cell_ptr->id() ) == _cell_id[ipf].end() ) continue;

    //filter out if not in energy range of sensitivity profile
    if ( ! within_energy_range( erg, ipf ) ) continue;

    //search for isotope in material
    std::shared_ptr< xs::neutronData > iso;
    double aden;
    {
      uint32_t i;
      bool     found = false;
      for ( i = 0 ; i < mat->size() ; ++i ) {
        iso = mat->xs(i);
        if ( iso->zaid() == _zaid[ipf] ) { found = true; aden = mat->aden(i); break; }
      }
      if ( ! found ) continue;
    }
    const auto kiso = iso->id();
    _iso_id_cache[ipf] = kiso;

    const auto ipe  = xs_buffer.urr_energy_index(   kiso );
    const auto rtp  = xs_buffer.urr_interp_factor(  kiso );
    const auto ilx  = xs_buffer.urr_left_xs_index(  kiso );
    const auto irx  = xs_buffer.urr_right_xs_index( kiso );

    const auto uxs  = xs_buffer.urr_table_xs( _mt[ipf], kiso );

    _pert_rxn_left_xs_cache[ipf]  = aden * uxs.first  * ( 1.0 - rtp );
    _pert_rxn_right_xs_cache[ipf] = aden * uxs.second * rtp;

    const auto ipl = ip + ipe    *(_nprob[ipf] + 1) + 1;
    const auto ipr = ip + (ipe+1)*(_nprob[ipf] + 1) + 1;

    _indirect_effect[ip]  -= wgt_times_distance * ( _pert_rxn_left_xs_cache[ipf] + _pert_rxn_right_xs_cache[ipf] );
    _indirect_effect[ipl] -= wgt_times_distance * _pert_rxn_left_xs_cache[ipf];
    _indirect_effect[ipr] -= wgt_times_distance * _pert_rxn_right_xs_cache[ipf];
    _indirect_effect[ipl + ilx + 1] -= wgt_times_distance * _pert_rxn_left_xs_cache[ipf];
    _indirect_effect[ipr + irx + 1] -= wgt_times_distance * _pert_rxn_right_xs_cache[ipf];
  }
}

void mc::dosUREstimator::score_indirect_collision( const double wgt, const double erg,
                                                   const std::shared_ptr<geom::cell> cell_ptr,
                                                   const xs::xsBuffer& xs_buffer ) {
  const auto mat = cell_ptr->material();
  const auto collided_zaid = xs_buffer.collided_zaid;
  const auto collided_mt   = xs_buffer.collided_mt;

  for ( auto ipf = 0 ; ipf < _nprofiles ; ++ipf ) {
    const auto ip = _parameter_id[ipf];

    //filter out if the cached value of the perturbed cross section is zero
    //usually implies where we are is not involved in the perturbation
    if ( _pert_rxn_left_xs_cache[ipf] == 0.0 && _pert_rxn_right_xs_cache[ipf] == 0.0 ) continue;

    //filter out not being in impacted cell
    //TODO: statement above means this filter could be removed
    if ( ! _cell_id[ipf].empty() && _cell_id[ipf].find( cell_ptr->id() ) == _cell_id[ipf].end() ) continue;

    //filter out reaction with isotope that does not match perturbation
    if ( _zaid[ipf] > 0 && _zaid[ipf] != collided_zaid ) continue;

    //filter out reaction that does not match perturbation
    if ( _mt[ipf] != xs::endf_mt_t::TOTAL && _mt[ipf] != collided_mt ) continue;

    //note: it is tempting to replace this with an implicit estimator scored every collision
    //      but unfortunately that is wrong because while the expectation of the indirect
    //      effect would be the same, what matters is the expectation of the PRODUCT of the
    //      indirect effect and the response score of the history, which would not be equal
    const auto kiso = _iso_id_cache[ipf];
    const auto ipe  = xs_buffer.urr_energy_index(   kiso );
    const auto rtp  = xs_buffer.urr_interp_factor(  kiso );
    const auto ilx  = xs_buffer.urr_left_xs_index(  kiso );
    const auto irx  = xs_buffer.urr_right_xs_index( kiso );

    const auto ipl = ip + ipe    *(_nprob[ipf] + 1) + 1;
    const auto ipr = ip + (ipe+1)*(_nprob[ipf] + 1) + 1;
//    const auto xs  = _pert_rxn_left_xs_cache[ipf] + _pert_rxn_right_xs_cache[ipf];
//    const auto fxl = _pert_rxn_left_xs_cache[ipf]  / xs;
//    const auto fxr = _pert_rxn_right_xs_cache[ipf] / xs;

    _indirect_effect[ip]  += wgt;
    _indirect_effect[ipl] += wgt * (1.0-rtp);
    _indirect_effect[ipr] += wgt * rtp;
    _indirect_effect[ipl + ilx + 1] += wgt * (1.0-rtp);
    _indirect_effect[ipr + irx + 1] += wgt * rtp;
  }
}

void mc::dosUREstimator::score_direct_response( const double wgt_times_distance, const double erg,
                                                const std::shared_ptr<geom::cell> cell_ptr,
                                                const xs::xsBuffer& xs_buffer ) {
  //grab material and filter out voids
  const auto mat = cell_ptr->material();
  if ( ! mat ) return;

  if ( _calc_mode == mc::mc_calc_t::K_EIGENVALUE ) {
    //score internal keff response first
    _estimators[_internal_keff_estimator_id]->score( wgt_times_distance, cell_ptr->id(), erg, cell_ptr, xs_buffer );
  }

  auto kit = 0;
  for ( auto eit = _estimators.begin() ; eit != _estimators.end() ; ++eit ) {
    const auto ir = _response_id[kit];
    ++kit;

    //direct effect is zero if no cross section assigned to the estimator
    if ( eit->second->mt() == xs::endf_mt_t::NONE ) continue;

    //filter out surface current estimators since direct effect is zero
    if ( eit->second->type() == mc::estimator_type_t::SURFACE_CURRENT ) continue;

    //filter out outside estimator response range
    if ( ! eit->second->within_energy_range( erg ) ) continue;
    const auto ire = eit->second->energy_bin( erg );

    const auto rzaid    = eit->second->zaid();
    const auto rmt      = eit->second->mt();
    for ( auto ipf = 0 ; ipf < _nprofiles ; ++ipf ) {
      const auto ip = _parameter_id[ipf];

      //filter out if estimator not within the perturbed cell or one scored by the estimator
      //TODO: verify this logic
      if ( ! _cell_id[ipf].empty() && _cell_id[ipf].find( cell_ptr->id() ) == _cell_id[ipf].end() ) continue;
      if ( ! eit->second->within_geom_set( cell_ptr->id() ) ) continue;

      //filter out if not in energy range of sensitivity profile
      if ( ! within_energy_range( erg, ipf ) ) continue;

      //filter out if the response isotope does not match perturbation if
      //a specific isotope was requested
      if ( rzaid != 0 && _zaid[ipf] != 0 && rzaid != _zaid[ipf] ) continue;

      //filter out if the response reaction does not match if not the total
      //note that nufission and fission need to be taken to be the same
      //if response is nufission and perturbation is fission we continue
      if ( rmt != xs::endf_mt_t::TOTAL && _mt[ipf] != xs::endf_mt_t::TOTAL && rmt != _mt[ipf] &&
        !( rmt == xs::endf_mt_t::NUFISSION && _mt[ipf] == xs::endf_mt_t::FISSION ) ) continue;

      //TODO: the total xs estimator multiplier is likely not handled correctly (also check base dosEstimator)
      //search for isotope in material and get atomic density
      double aden;
      {
        std::shared_ptr< xs::neutronData > iso;
        uint32_t i;
        bool     found = false;
        for ( i = 0 ; i < mat->size() ; ++i ) {
          iso = mat->xs(i);
          if ( iso->zaid() == _zaid[ipf] ) { found = true; aden = mat->aden(i); break; }
        }
        if ( ! found ) continue;
      }

      const auto kiso = _iso_id_cache[ipf];
      const auto ipe  = xs_buffer.urr_energy_index(   kiso );
      const auto rtp  = xs_buffer.urr_interp_factor(  kiso );
      const auto ilx  = xs_buffer.urr_left_xs_index(  kiso );
      const auto irx  = xs_buffer.urr_right_xs_index( kiso );

      const auto uxs  = xs_buffer.urr_table_xs( rmt, kiso );
      const auto xsl = aden * uxs.first  * ( 1.0 - rtp );
      const auto xsr = aden * uxs.second * rtp;

      const auto ipl = ip + ipe    *(_nprob[ipf] + 1) + 1;
      const auto ipr = ip + (ipe+1)*(_nprob[ipf] + 1) + 1;

      const auto sl = wgt_times_distance * xsl;
      const auto sr = wgt_times_distance * xsr;

      _direct_effect[ir][ip]  += sl + sr;
      _direct_effect[ir][ipl] += sl;
      _direct_effect[ir][ipr] += sr;
      _direct_effect[ir][ipl + ilx + 1] += sl;
      _direct_effect[ir][ipr + irx + 1] += sr;
      if ( ire > 0 ) {
        _direct_effect[ir+ire][ip]  += sl + sr;
        _direct_effect[ir+ire][ipl] += sl;
        _direct_effect[ir+ire][ipr] += sr;
        _direct_effect[ir+ire][ipl + ilx + 1] += sl;
        _direct_effect[ir+ire][ipr + irx + 1] += sr;
      }
    }
  }
}

//bank perturbation for collection of fission source points
//called once when any fission neutrons are banked
uint32_t mc::dosUREstimator::insert_fission_source( const double wgt, const double erg,
                                                    const std::shared_ptr<geom::cell> cell_ptr,
                                                    const uint32_t num_neutrons,
                                                    const xs::xsBuffer& xs_buffer ) {
  assert( _calc_mode == mc::mc_calc_t::K_EIGENVALUE );

  //compute direct effect of perturbation on making fission src particle with collision estimator
  const auto mat = cell_ptr->material();

  for ( auto ipf = 0 ; ipf < _nprofiles ; ++ipf ) {
    const auto ip  = _parameter_id[ipf];

    double direct_effect_nufission       = 0.0;
    double direct_effect_nufission_left  = 0.0;
    double direct_effect_nufission_right = 0.0;
    if ( _pert_rxn_left_xs_cache[ipf] != 0.0 || _pert_rxn_right_xs_cache[ipf] != 0.0 ) {
      //if the reaction xs > 0 then isotope in the same material
      double total_xs   = 0.0;
      double fission_xs = 0.0;
      double df         = 0.0;
      for ( auto i = 0 ; i < mat->size() ; ++i ) {
        const auto iso = mat->xs(i);
        const auto k   = iso->id();
        total_xs   += mat->aden(i) * xs_buffer.total_xs[k];
        fission_xs += mat->aden(i) * xs_buffer.fission_xs[k];
      }
      df -= 1.0/total_xs;  //TODO: handle not a cross section case...
      //total in this context can be thought of as a density perturbation for either
      //the entire material or a single isotope
      if ( _mt[ipf] == xs::endf_mt_t::FISSION || _mt[ipf] == xs::endf_mt_t::TOTAL ) {
        df += 1.0/fission_xs;
      }
      direct_effect_nufission_left  = wgt * df * _pert_rxn_left_xs_cache[ipf];
      direct_effect_nufission_right = wgt * df * _pert_rxn_right_xs_cache[ipf];
      direct_effect_nufission       = direct_effect_nufission_left + direct_effect_nufission_right;
    }
    const auto pert_wgt = _indirect_effect[ip] + direct_effect_nufission;

    //store perturbed weight and accrue fission source normalization for total
    _fission_bank[_fbank_index][ip]  = pert_wgt;
    _total_perturbed_weight[ip]     += pert_wgt * num_neutrons;

    //loop over energy bins
    const auto kiso = _iso_id_cache[ipf];
    const auto ipe  = xs_buffer.urr_energy_index(   kiso );
    const auto ilx  = xs_buffer.urr_left_xs_index(  kiso );
    const auto irx  = xs_buffer.urr_right_xs_index( kiso );
    auto k = 1;
    for ( auto ie = 0 ; ie < _nebins[ipf] ; ++ie ) {
      auto pert_wgt = _indirect_effect[ip+k];
      pert_wgt += ie == ipe   ? direct_effect_nufission_left  : 0.0;
      pert_wgt += ie == ipe+1 ? direct_effect_nufission_right : 0.0;

      _fission_bank[_fbank_index][ip+k]  = pert_wgt;
      _total_perturbed_weight[ip+k]     += pert_wgt * num_neutrons;
      ++k;

      //loop over probability levels
      for ( auto ix = 0 ; ix < _nprob[ipf] ; ++ix ) {
        auto pert_wgt = _indirect_effect[ip+k];
        pert_wgt += ie == ipe   && ix == ilx ? direct_effect_nufission_left  : 0.0;
        pert_wgt += ie == ipe+1 && ix == irx ? direct_effect_nufission_right : 0.0;

        _fission_bank[_fbank_index][ip+k]  = pert_wgt;
        _total_perturbed_weight[ip+k]     += pert_wgt * num_neutrons;
        ++k;
      }
    }
  }
  return _fbank_index++;
}

void mc::dosUREstimator::end_history() {
  auto kr = 0;
  for ( auto eit = _estimators.begin() ; eit != _estimators.end() ; ++eit ) {
    const auto ir  = _response_id[kr];
    const auto eid = eit->second->id();
    const auto nre = eit->second->nebins() >= 2 ? eit->second->nebins() + 1 : 1;

    //loop over all responses contained in estimator
    for ( auto ire = 0 ; ire < nre ; ++ire ) {
      //loop over all parameters
      for ( auto ip = 0 ; ip < _nparameters ; ++ip ) {
        const auto score = _direct_effect[ir+ire][ip] + _indirect_effect[ip] * eit->second->current_score(ire);
        _sens_mean[eid][ire][ip] += score;
        _sens_stdv[eid][ire][ip] += score * score;
      }
    }
    eit->second->end_history();
    ++kr;
  }
}

void mc::dosUREstimator::end_batch( const uint32_t total_neutrons_banked ) {
  if ( _calc_mode == mc::mc_calc_t::K_EIGENVALUE ) {
    for ( auto ip = 0 ; ip < _nparameters ; ++ip ) {
      _fission_source_renorm[ip]  = _total_perturbed_weight[ip] / total_neutrons_banked;
      _total_perturbed_weight[ip] = 0.0;
    }
    _fbank_index = 0;
    std::swap( _fission_bank, _source_bank );
  }
}

void mc::dosUREstimator::write_sensitivity( const double normalization ) const {
  if ( _calc_mode == mc::mc_calc_t::K_EIGENVALUE ) {
    //write out keff sensitivity first
    const auto keff_mean = _estimators.find(_internal_keff_estimator_id)->second->mean( normalization, 0 );
    const auto keff_stdv = _estimators.find(_internal_keff_estimator_id)->second->stdv( normalization, 0 );

    std::cout << '\n' << " keff = ";
    std::cout << std::fixed << std::setprecision(8) << keff_mean << " +/- " << keff_stdv << '\n';
    for ( auto ipf = 0 ; ipf < _nprofiles ; ++ipf ) {
      const auto ip = _parameter_id[ipf];

      std::cout << " sensitivity for " << _name[ipf] << " = " << '\n';
      for ( auto ipe = 0 ; ipe < _nebins[ipf] ; ++ipe ) {
        std::cout << std::scientific << std::setprecision(3) << _egrid[ipf][ipe] << ' ';

        double sum_lev = 0.0;  //TODO: remove this debugging variable

        const auto ke = ip + ipe*(_nprob[ipf] + 1) + 1;
        for ( auto ipx = 0 ; ipx < _nprob[ipf] ; ++ipx ) {
          const auto ksen_mean = _sens_mean.find(_internal_keff_estimator_id)->second[0][ke+ipx+1] / normalization;
          std::cout << std::scientific << std::setprecision(3) << std::setw(10) << ksen_mean / keff_mean << ' ';
          sum_lev += ksen_mean / keff_mean;
        }
        //print out sum over probability levels at end of row
        const auto ksen_mean = _sens_mean.find(_internal_keff_estimator_id)->second[0][ke] / normalization;
        std::cout << std::scientific << std::setprecision(3) << std::setw(10) << ksen_mean / keff_mean << '\n';
        assert( std::fabs( sum_lev - ksen_mean / keff_mean ) <= 1.0e-6 );
      }
      const auto ksen_mean = _sens_mean.find(_internal_keff_estimator_id)->second[0][ip] / normalization;
      //TODO: compute better uncertainty estimate (should include correlations)
      const auto ksen2     = _sens_stdv.find(_internal_keff_estimator_id)->second[0][ip] / normalization;
      const auto ksen_stdv = std::sqrt( ( ksen2 - ksen_mean*ksen_mean )/normalization );

      std::cout << " Total = " << std::scientific << std::setprecision(6) << ksen_mean / keff_mean
                << " +/- " << ksen_stdv / keff_mean << '\n';
      std::cout << std::fixed;
    }
  }

  //find maximum string length
  size_t max_strlen = 0;
  for ( auto rrit : _rxn_ratios ) {
    const auto numerator    = rrit.first;
    const auto nest_name = _estimators.find( numerator->id() )->second->name();

    const auto denominator  = rrit.second;
    if ( denominator ) {
      const auto dest_name = _estimators.find( denominator->id() )->second->name();
      std::string str = nest_name + " / " + dest_name;
      max_strlen = std::max( max_strlen, str.length() );
    }
    else {
      max_strlen = std::max( max_strlen, nest_name.length() );
    }
  }

  //then write out all user-defined reaction rate ratio sensitivities
  for ( auto rrit : _rxn_ratios ) {
    const auto numerator    = rrit.first;
    const auto denominator  = rrit.second;

    //number of bins in numerator and denominator
    //denominator must either have one bin or an equal number to numerator
    const auto nbins = rrit.first->nebins()  >= 2 ? rrit.first->nebins()  + 1 : 1;
    const auto dbins = denominator ? ( rrit.second->nebins() >= 2 ? rrit.second->nebins() + 1 : 1 ) : 1;
    assert( dbins == 1 || dbins == nbins );

    std::string rr_str;
    const auto& nest = _estimators.find( numerator->id() )->second;
    const auto nest_name = nest->name();
    if ( denominator ) {
      const auto& dest = _estimators.find( denominator->id() )->second;
      const auto dest_name = dest->name();

      rr_str = nest_name + " / " + dest_name;
    }
    else {
      rr_str = nest_name;
    }

    std::cout << '\n' << std::setw(max_strlen) << rr_str;
    for ( auto inbin = 0 ; inbin < nbins ; ++inbin ) {
      const auto nest_mean = nest->mean( normalization, inbin );
      const auto nest_stdv = nest->stdv( normalization, inbin );

      const auto idbin = dbins == 1 ? 0 : inbin;
      const auto dest_mean = denominator ? _estimators.find( denominator->id() )->second->mean( normalization, idbin ) : 1.0;
      const auto dest_stdv = denominator ? _estimators.find( denominator->id() )->second->stdv( normalization, idbin ) : 0.0;

      const auto rr_mean = nest_mean / dest_mean;
      const auto rr_stdv = rr_mean * std::sqrt( std::pow( nest_stdv/nest_mean, 2 )
                                              + std::pow( dest_stdv/dest_mean, 2 ) );
      std::cout << std::scientific;
      if ( inbin == 0 ) {
        std::cout << " = " << std::setprecision(6) << std::setw(12)
                  << rr_mean << " +/- " << rr_stdv << '\n';
      }
      else {
        std::cout << '\n' << rr_str << " from " << std::setprecision(4)
                  << nest->egrid_value(inbin-1) << " to " << nest->egrid_value(inbin) << " = "
                  << std::setprecision(6) << std::setw(12)
                  << rr_mean << " +/- " << rr_stdv << '\n';
      }
      for ( auto ipf = 0 ; ipf < _nprofiles ; ++ipf ) {
        const auto ip = _parameter_id[ipf];

        std::cout << " sensitivity for " << _name[ipf] << " = " << '\n';
        for ( auto ipe = 0 ; ipe < _nebins[ipf] ; ++ipe ) {
          std::cout << std::scientific << std::setprecision(3) << _egrid[ipf][ipe] << ' ';

          double sum_lev = 0.0;  //TODO: remove this debugging variable

          const auto ke = ip + ipe*(_nprob[ipf] + 1) + 1;
          for ( auto ipx = 0 ; ipx < _nprob[ipf] ; ++ipx ) {
            const auto nsen_mean = _sens_mean.find( numerator->id()   )->second[inbin][ke+ipx+1] / normalization;
            const auto dsen_mean = denominator ? _sens_mean.find( denominator->id() )->second[idbin][ke+ipx+1] / normalization : 0.0;
            const auto rr_sen_mean = nsen_mean / nest_mean - dsen_mean / dest_mean;

            std::cout << std::scientific << std::setprecision(3) << std::setw(10) << rr_sen_mean << ' ';
            sum_lev += rr_sen_mean;
          }
          //print out sum over probability levels at end of row
          const auto nsen_mean = _sens_mean.find( numerator->id()   )->second[inbin][ke] / normalization;
          const auto dsen_mean = denominator ? _sens_mean.find( denominator->id() )->second[idbin][ke] / normalization : 0.0;
          const auto rr_sen_mean = nsen_mean / nest_mean - dsen_mean / dest_mean;

          std::cout << std::scientific << std::setprecision(3) << std::setw(10) << rr_sen_mean << '\n';
          assert( std::fabs( sum_lev - rr_sen_mean ) <= 1.0e-6 );
        }
        //print out total with uncertainty
        const auto nsen_mean = _sens_mean.find( numerator->id() )->second[inbin][ip] / normalization;
        const auto nsen2     = _sens_stdv.find( numerator->id() )->second[inbin][ip] / normalization;
        const auto nsen_var  = ( nsen2 - nsen_mean*nsen_mean ) / normalization;
        const auto nsen_stdv = std::sqrt( nsen_var );
        const auto nvar = std::pow( nsen_mean / nest_mean, 2 )
                            * ( std::pow( nsen_stdv/nsen_mean, 2 ) + std::pow( nest_stdv/nest_mean, 2 )  );

        double dsen_mean = 0.0;
        double dvar      = 0.0;
        if ( denominator ) {
          dsen_mean = _sens_mean.find( denominator->id() )->second[idbin][ip] / normalization;

          const auto dsen2     = _sens_stdv.find( denominator->id() )->second[idbin][ip] / normalization;
          const auto dsen_var  = ( dsen2 - dsen_mean*dsen_mean ) / normalization;
          const auto dsen_stdv = std::sqrt( dsen_var );
          dvar = std::pow( dsen_mean / dest_mean, 2 )
               * ( std::pow( dsen_stdv/dsen_mean, 2 ) + std::pow( dest_stdv/dest_mean, 2 )  );
        }

        const auto rr_sen_mean = nsen_mean / nest_mean - dsen_mean / dest_mean;
        const auto rr_sen_stdv = std::sqrt( nvar + dvar );

        std::cout << " Total = " << std::scientific << std::setprecision(6) << rr_sen_mean
                  << " +/- " << rr_sen_stdv << '\n';
/*
        const auto ip = _parameter_id[ipf];

        const auto nsen_mean = _sens_mean.find( numerator->id() )->second[inbin][ip] / normalization;
        const auto nsen2     = _sens_stdv.find( numerator->id() )->second[inbin][ip] / normalization;
        const auto nsen_var  = ( nsen2 - nsen_mean*nsen_mean ) / normalization;
        const auto nsen_stdv = std::sqrt( nsen_var );

        const auto dsen_mean = _sens_mean.find( denominator->id() )->second[idbin][ip] / normalization;
        const auto dsen2     = _sens_stdv.find( denominator->id() )->second[idbin][ip] / normalization;
        const auto dsen_var  = ( dsen2 - dsen_mean*dsen_mean ) / normalization;
        const auto dsen_stdv = std::sqrt( dsen_var );

        const auto nvar = std::pow( nsen_mean / nest_mean, 2 )
                            * ( std::pow( nsen_stdv/nsen_mean, 2 ) + std::pow( nest_stdv/nest_mean, 2 )  );
        const auto dvar = std::pow( dsen_mean / dest_mean, 2 )
                            * ( std::pow( dsen_stdv/dsen_mean, 2 ) + std::pow( dest_stdv/dest_mean, 2 )  );

        const auto rr_sen_mean = nsen_mean / nest_mean - dsen_mean / dest_mean;
        const auto rr_sen_stdv = std::sqrt( nvar + dvar );

        std::cout << " sensitivity for " << _name[ipf] << " = ";
        if ( _nebins[ipf] == 0 ) {
          //single parameter bin
          std::cout << std::setprecision(8) << rr_sen_mean << " +/- " << rr_sen_stdv << '\n';
        }
        else {
          //multiple parameter bins
          std::cout << '\n';
          for ( auto ipe = 1 ; ipe < _egrid[ipf].size() ; ++ipe ) {
            const auto nsen_mean = _sens_mean.find( numerator->id() )->second[inbin][ip+ipe] / normalization;
            const auto nsen2     = _sens_stdv.find( numerator->id() )->second[inbin][ip+ipe] / normalization;
            const auto nsen_var  = ( nsen2 - nsen_mean*nsen_mean ) / normalization;
            const auto nsen_stdv = std::sqrt( nsen_var );

            const auto dsen_mean = _sens_mean.find( denominator->id() )->second[idbin][ip+ipe] / normalization;
            const auto dsen2     = _sens_stdv.find( denominator->id() )->second[idbin][ip+ipe] / normalization;
            const auto dsen_var  = ( dsen2 - dsen_mean*dsen_mean ) / normalization;
            const auto dsen_stdv = std::sqrt( dsen_var );

            const auto nvar = std::pow( nsen_mean / nest_mean, 2 )
                                * ( std::pow( nsen_stdv/nsen_mean, 2 ) + std::pow( nest_stdv/nest_mean, 2 )  );
            const auto dvar = std::pow( dsen_mean / dest_mean, 2 )
                                * ( std::pow( dsen_stdv/dsen_mean, 2 ) + std::pow( dest_stdv/dest_mean, 2 )  );

            const auto rr_sen_mean = nsen_mean / nest_mean - dsen_mean / dest_mean;
            const auto rr_sen_stdv = std::sqrt( nvar + dvar );

            std::cout << std::scientific << std::setprecision(4) << std::setw(10)
                      << _egrid[ipf][ipe-1] << " to " << _egrid[ipf][ipe] << " = "
                      << std::setprecision(6) << std::setw(12)
                      << rr_sen_mean << " +/- " << rr_sen_stdv << '\n';
          }
          std::cout << " Total =  " << std::setprecision(8)
                    << rr_sen_mean << " +/- " << rr_sen_stdv << '\n';
        }
*/
      }
    }
  }
}
