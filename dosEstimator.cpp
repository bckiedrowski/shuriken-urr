#include <cmath>
#include <iomanip>
#include <iostream>

#include "binarySearch.hpp"

#include "dosEstimator.hpp"

mc::dosEstimator::dosEstimator( const std::vector< dos_profile_t > profiles,
                                const std::optional< rxn_ratio_vector > rxn_ratios,
                                const mc::mc_calc_t calc_mode ) {
  _calc_mode = calc_mode;
  _nprofiles = profiles.size();

  //number of parameters is the sum of the size of each of the profiles
  _nparameters = 0;
  for ( auto pr : profiles ) {
    if ( pr.egrid ) {
      //number of parameters in this profile is the number of energy bins plus one for the total
      //except when a single energy bin is given, then there is only one parameter
      const auto& egrid = pr.egrid.value();
      assert( egrid.size() >= 2 );  //need at least two grid points for a single bin
      if ( egrid.size() > 2 ) {
        //number of parameters is egrid.size() - 1 for energy grid and + 1 for total bin
        _nparameters += egrid.size();
      }
      else {
        //single energy bin is special case for one response since there is no point to
        //compute for both a single bin and a sum total over just that bin
        _nparameters++;
      }
    }
    else {
      //if no energy grid given then one parameter in this profile
      _nparameters++;
    }
  }

  _name.resize( _nprofiles );
  _cell_id.resize( _nprofiles );
  _zaid.resize( _nprofiles );
  _mt.resize( _nprofiles );
  _nebins.resize( _nprofiles );
  _egrid.resize( _nprofiles );
  _parameter_id.resize( _nprofiles );
  _within_xs_buffer.resize( _nprofiles );

  _pert_rxn_xs_cache.resize( _nprofiles, 0.0 );
  _erg_offset_cache.resize( _nprofiles, 0 );

  auto ip = 0;   //ip = parameter index, ipf = profile index
  for ( auto ipf = 0 ; ipf < _nprofiles ; ++ipf ) {
    _name[ipf]    = profiles[ipf].name;
    _cell_id[ipf] = profiles[ipf].cell_id.value_or( std::set<uint32_t>() );
    _zaid[ipf]    = profiles[ipf].zaid.value_or( 0 );
    _mt[ipf]      = profiles[ipf].mt.value_or( xs::endf_mt_t::TOTAL );

    if ( profiles[ipf].egrid ) {
      const auto& egrid = profiles[ipf].egrid.value();
      _egrid[ipf]  = egrid;
      _nebins[ipf] = egrid.size() > 2 ? egrid.size()-1 : 1;
    }
    else {
      _nebins[ipf] = 0;  //denotes over all energies
    }

    _parameter_id[ipf] = ip;
    //TODO: verify this
    ip += std::max( _nebins[ipf] + 1, 1u );

    //if zaid = 0, this is a density perturbation and it only makes sense to specify
    //the total cross section
    const auto& mt = _mt[ipf];
    if ( _zaid[ipf] == 0 ) assert( mt == xs::endf_mt_t::TOTAL );

    _within_xs_buffer[ipf] =
      ( mt == xs::endf_mt_t::TOTAL   || mt == xs::endf_mt_t::ELASTIC ||
        mt == xs::endf_mt_t::FISSION || mt == xs::endf_mt_t::Z_GAMMA    );
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

bool mc::dosEstimator::within_energy_range( const double erg, const uint32_t ipf ) const {
  if ( _nebins[ipf] > 0 ) {
    if ( erg < _egrid[ipf][0] || erg >= _egrid[ipf][_nebins[ipf]] ) return false;
  }
  return true;
}

void mc::dosEstimator::resize_source_bank( const uint32_t bank_size ) {
  assert( _calc_mode == mc::mc_calc_t::K_EIGENVALUE );

  _fission_source_renorm.resize( _nparameters, 0.0 );
  _total_perturbed_weight.resize( _nparameters, 0.0 );

  _fission_bank.resize( bank_size, std::vector< double >( _nparameters, 0.0 ) );
  _source_bank.resize(  bank_size, std::vector< double >( _nparameters, 0.0 ) );
}

void mc::dosEstimator::start_history( const uint32_t bank_index ) {
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

void mc::dosEstimator::score_indirect_stream( const double wgt_times_distance, const double erg,
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

    //filter out not being in impacted cell
    if ( ! _cell_id[ipf].empty() && _cell_id[ipf].find( cell_ptr->id() ) == _cell_id[ipf].end() ) continue;

    //filter out if not in energy range of sensitivity profile
    if ( ! within_energy_range( erg, ipf ) ) continue;

    //compute and store energy index offset (remains zero if single bin)
    uint32_t ipe = 0;
    if ( _nebins[ipf] >= 2 ) {
      if ( _nebins[ipf] < 10 ) {
        //linear search for small energy grid
        for ( ipe = 1 ; ipe < _nebins[ipf] ; ++ipe ) {
          if ( _egrid[ipf][ipe] > erg ) break;
        }
      }
      else {
        ipe = 1 + util::binary_search( erg, _egrid[ipf] );
      }
      _erg_offset_cache[ipf] = ipe;
    }

    //compute and store perturbed cross section
    if ( _zaid[ipf] == 0 ) {
      //TODO: this only makes sense for the total cross section
      //need cross section for entire material since all isotopes in perturbation
      for ( auto i = 0 ; i < mat->size() ; ++i ) {
        const auto iso = mat->xs(i);
        const auto k   = iso->id();
        _pert_rxn_xs_cache[ipf] += mat->aden(i) *
          ( _within_xs_buffer[ipf] ? xs_buffer( _mt[ipf], k ) : iso->find_xs( _mt[ipf], erg ) );
      }
    }
    else {
      //only for the specific isotope requested
      //need to check if the isotope is in the current material
      std::shared_ptr< xs::neutronData > iso;
      uint32_t i;
      bool     found = false;
      for ( i = 0 ; i < mat->size() ; ++i ) {
        iso = mat->xs(i);
        if ( iso->zaid() == _zaid[ipf] ) { found = true; break; }
      }
      if ( found ) {
        const auto k  = iso->id();
        _pert_rxn_xs_cache[ipf] += mat->aden(i) *
          ( _within_xs_buffer[ipf] ? xs_buffer( _mt[ipf], k ) : iso->find_xs( _mt[ipf], erg ) );
      }
    }
    _indirect_effect[ip] -= wgt_times_distance * _pert_rxn_xs_cache[ipf];
    if ( ipe > 0 ) _indirect_effect[ip+ipe] -= wgt_times_distance * _pert_rxn_xs_cache[ipf];
  }
}

void mc::dosEstimator::score_indirect_collision( const double wgt, const double erg,
                                                 const std::shared_ptr<geom::cell> cell_ptr,
                                                 const xs::xsBuffer& xs_buffer ) {
  const auto mat = cell_ptr->material();
  const auto collided_zaid = xs_buffer.collided_zaid;
  const auto collided_mt   = xs_buffer.collided_mt;

  for ( auto ipf = 0 ; ipf < _nprofiles ; ++ipf ) {
    const auto ip = _parameter_id[ipf];

    //filter out if the cached value of the perturbed cross section is zero
    //usually implies where we are is not involved in the perturbation
    if ( _pert_rxn_xs_cache[ipf] == 0.0 ) continue;

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
    _indirect_effect[ip] += wgt;

    const auto ipe = _erg_offset_cache[ipf];
    if ( ipe > 0 ) _indirect_effect[ip+ipe] += wgt;
  }
}

void mc::dosEstimator::score_direct_response( const double wgt_times_distance, const double erg,
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
      const auto ipe = _erg_offset_cache[ipf];

      //filter out if the response isotope does not match perturbation if
      //a specific isotope was requested
      if ( rzaid != 0 && _zaid[ipf] != 0 && rzaid != _zaid[ipf] ) continue;

      //filter out if the response reaction does not match if not the total
      //note that nufission and fission need to be taken to be the same
      //if response is nufission and perturbation is fission we continue
      if ( rmt != xs::endf_mt_t::TOTAL && _mt[ipf] != xs::endf_mt_t::TOTAL && rmt != _mt[ipf] &&
        !( rmt == xs::endf_mt_t::NUFISSION && _mt[ipf] == xs::endf_mt_t::FISSION ) ) continue;

      if ( _zaid[ipf] == 0 ) {
        //density perturbation over the all isotopes
        if ( rzaid == 0 ) {
          //response is for all isotopes in the current cell for response reaction
          auto rxn_xs = 0.0;
          for ( auto i = 0 ; i < mat->size() ; ++i ) {
            const auto iso = mat->xs(i);
            const auto k   = iso->id();
            rxn_xs += mat->aden(i) *
              ( _within_xs_buffer[ipf] ? xs_buffer( rmt, k ) : iso->find_xs( rmt, erg ) );
          }
          const auto s = wgt_times_distance * rxn_xs;
          _direct_effect[ir][ip] += s;
          if ( ipe > 0 ) _direct_effect[ir][ip+ipe] += s;
          if ( ire > 0 ) _direct_effect[ir+ire][ip] += s;
          if ( ipe > 0 && ire > 0 ) _direct_effect[ir+ire][ip+ipe] += s;
        }
        else {
          //response for a single isotope (match confirmed above) [TODO: verify]
          std::shared_ptr< xs::neutronData > iso;
          uint32_t i;
          bool     found_match = false;
          for ( i = 0 ; i < mat->size() ; ++i ) {
            iso = mat->xs(i);
            if ( iso->zaid() == rzaid ) { found_match = true; break; }
          }
          if ( found_match ) {
            const auto k  = iso->id();
            const auto s  = wgt_times_distance * mat->aden(i) *
              ( _within_xs_buffer[ipf] ? xs_buffer( rmt, k ) : iso->find_xs( rmt, erg ) );
            _direct_effect[ir][ip] += s;
            if ( ipe > 0 ) _direct_effect[ir][ip+ipe] += s;
            if ( ire > 0 ) _direct_effect[ir+ire][ip] += s;
            if ( ipe > 0 && ire > 0 ) _direct_effect[ir+ire][ip+ipe] += s;
          }
        }
      }
      else {
        //only a specific isotope is being perturbed
        //need to check if the isotope is in the current material
        std::shared_ptr< xs::neutronData > iso;
        uint32_t i;
        bool     found_match = false;
        for ( i = 0 ; i < mat->size() ; ++i ) {
          iso = mat->xs(i);
          if ( iso->zaid() == _zaid[ipf] ) { found_match = true; break; }
        }
        if ( found_match ) {
          const auto k  = iso->id();
          const auto s  = wgt_times_distance * mat->aden(i) *
            ( _within_xs_buffer[ipf] ? xs_buffer( rmt, k ) : iso->find_xs( rmt, erg ) );
          _direct_effect[ir][ip] += s;
          if ( ipe > 0 ) _direct_effect[ir][ip+ipe] += s;
          if ( ire > 0 ) _direct_effect[ir+ire][ip] += s;
          if ( ipe > 0 && ire > 0 ) _direct_effect[ir+ire][ip+ipe] += s;
        }
      }
    }
  }
}

//bank perturbation for collection of fission source points
//called once when any fission neutrons are banked
uint32_t mc::dosEstimator::insert_fission_source( const double wgt, const double erg,
                                                  const std::shared_ptr<geom::cell> cell_ptr,
                                                  const uint32_t num_neutrons,
                                                  const xs::xsBuffer& xs_buffer ) {
  assert( _calc_mode == mc::mc_calc_t::K_EIGENVALUE );

  //compute direct effect of perturbation on making fission src particle with collision estimator
  const auto mat = cell_ptr->material();

  for ( auto ipf = 0 ; ipf < _nprofiles ; ++ipf ) {
    const auto ip  = _parameter_id[ipf];
    const auto ipe = _erg_offset_cache[ipf];

    double direct_effect_nufission = 0.0;
    if ( _pert_rxn_xs_cache[ipf] > 0.0 ) {
      //if the reaction xs > 0 then isotope in the same material
      //TODO: handle nubar
      //TODO: handle chances of fission...
      //TODO: duplicate code with above, need to clean up
      double total_xs   = 0.0;
      double fission_xs = 0.0;
      for ( auto i = 0 ; i < mat->size() ; ++i ) {
        const auto iso = mat->xs(i);
        const auto k   = iso->id();
        total_xs   += mat->aden(i) * xs_buffer.total_xs[k];
        fission_xs += mat->aden(i) * xs_buffer.fission_xs[k];
      }
      direct_effect_nufission -= 1.0/total_xs;  //TODO: handle not a cross section case...
      //total in this context can be thought of as a density perturbation for either
      //the entire material or a single isotope
      if ( _mt[ipf] == xs::endf_mt_t::FISSION || _mt[ipf] == xs::endf_mt_t::TOTAL ) {
        direct_effect_nufission += 1.0/fission_xs;
      }
      direct_effect_nufission *= wgt * _pert_rxn_xs_cache[ipf];
    }
    //TODO: verify that we do not need to multiply indirect effect by nu*sigma_f/sigma_t
    //      might make sense considering that we bank multiple neutrons
    const auto pert_wgt = _indirect_effect[ip] + direct_effect_nufission;

    //store perturbed weight and accrue fission source normalization for total
    _fission_bank[_fbank_index][ip]  = pert_wgt;
    _total_perturbed_weight[ip]     += pert_wgt * num_neutrons;

    //if multiple energy bins are provided, then bank for each energy bin
    if ( _nebins[ipf] > 1 ) {
      for ( auto k = 1 ; k <= _nebins[ipf] ; ++k ) {
        //direct effect only applies to energy bin ipe for the currrent particle
        const auto pert_wgt = _indirect_effect[ip+k] + ( k == ipe ? direct_effect_nufission : 0.0 );
        _fission_bank[_fbank_index][ip+k]  = pert_wgt;
        _total_perturbed_weight[ip+k]     += pert_wgt * num_neutrons;
      }
    }
  }
  return _fbank_index++;
}

void mc::dosEstimator::end_history() {
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

void mc::dosEstimator::end_batch( const uint32_t total_neutrons_banked ) {
  if ( _calc_mode == mc::mc_calc_t::K_EIGENVALUE ) {
    for ( auto ip = 0 ; ip < _nparameters ; ++ip ) {
      _fission_source_renorm[ip]  = _total_perturbed_weight[ip] / total_neutrons_banked;
      _total_perturbed_weight[ip] = 0.0;
    }
    _fbank_index = 0;
    std::swap( _fission_bank, _source_bank );
  }
}

void mc::dosEstimator::write_sensitivity( const double normalization ) const {
  if ( _calc_mode == mc::mc_calc_t::K_EIGENVALUE ) {
    //write out keff sensitivity first
    const auto keff_mean = _estimators.find(_internal_keff_estimator_id)->second->mean( normalization, 0 );
    const auto keff_stdv = _estimators.find(_internal_keff_estimator_id)->second->stdv( normalization, 0 );

    std::cout << '\n' << " keff = ";
    std::cout << std::fixed << std::setprecision(8) << keff_mean << " +/- " << keff_stdv << '\n';
    for ( auto ipf = 0 ; ipf < _nprofiles ; ++ipf ) {
      const auto ip = _parameter_id[ipf];

      std::cout << " sensitivity for " << _name[ipf] << " = ";
      const auto ksen_mean = _sens_mean.find(_internal_keff_estimator_id)->second[0][ip] / normalization;

      //TODO: compute better uncertainty estimate (should include correlations)
      const auto ksen2     = _sens_stdv.find(_internal_keff_estimator_id)->second[0][ip] / normalization;
      const auto ksen_stdv = std::sqrt( ( ksen2 - ksen_mean*ksen_mean )/normalization );
      if ( _nebins[ipf] == 0 ) {
        std::cout << std::setprecision(8) << ksen_mean / keff_mean << " +/- " << ksen_stdv / keff_mean << '\n';
      }
      else {
        std::cout << '\n';
        for ( auto ipe = 1 ; ipe < _egrid[ipf].size() ; ++ipe ) {
          const auto ksen_mean = _sens_mean.find(_internal_keff_estimator_id)->second[0][ip+ipe] / normalization;

          //TODO: compute better uncertainty estimate (should include correlations)
          const auto ksen2     = _sens_stdv.find(_internal_keff_estimator_id)->second[0][ip+ipe] / normalization;
          const auto ksen_stdv = std::sqrt( ( ksen2 - ksen_mean*ksen_mean )/normalization );

          std::cout << std::scientific << std::setprecision(4) << std::setw(10)
                    << _egrid[ipf][ipe-1] << " to " << _egrid[ipf][ipe] << " = "
                    << std::setprecision(6) << std::setw(12)
                    << ksen_mean / keff_mean << " +/- " << ksen_stdv / keff_mean << '\n';
        }
        std::cout << " Total =  " << std::setprecision(8)
                  << ksen_mean / keff_mean << " +/- " << ksen_stdv / keff_mean << '\n';
      }
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
            const auto nvar = std::pow( nsen_mean / nest_mean, 2 )
                                * ( std::pow( nsen_stdv/nsen_mean, 2 ) + std::pow( nest_stdv/nest_mean, 2 )  );

            double dsen_mean = 0.0;
            double dvar      = 0.0;
            if ( denominator ) {
              dsen_mean = _sens_mean.find( denominator->id() )->second[idbin][ip+ipe] / normalization;

              const auto dsen2     = _sens_stdv.find( denominator->id() )->second[idbin][ip+ipe] / normalization;
              const auto dsen_var  = ( dsen2 - dsen_mean*dsen_mean ) / normalization;
              const auto dsen_stdv = std::sqrt( dsen_var );
              dvar = std::pow( dsen_mean / dest_mean, 2 )
                   * ( std::pow( dsen_stdv/dsen_mean, 2 ) + std::pow( dest_stdv/dest_mean, 2 )  );
            }

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
      }
    }
  }
}
