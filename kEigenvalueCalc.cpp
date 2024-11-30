#include <cmath>
#include <iomanip>
#include <iostream>

#include "constants.hpp"

#include "kEigenvalueCalc.hpp"

mc::kEigenvalueCalc::kEigenvalueCalc( const size_t batch_size, const size_t num_inactive, const size_t num_batches,
                                      const util::point pos, const std::vector< std::shared_ptr<geom::cell> >& cells,
                                      const std::vector< std::shared_ptr< mc::estimator > >& estimators,
                                      const std::optional< std::shared_ptr< mc::dosEstimator > >   dos_estimators,
                                      const std::optional< std::shared_ptr< mc::dosUREstimator > > dos_ur_estimators )
  : _num_inactive(num_inactive), _batch_size_initial(batch_size),
    monteCarloCalc( batch_size, num_batches, estimators, dos_estimators, dos_ur_estimators ) {

   _source_bank.resize(  1.5 * _batch_size_initial );
   _fission_bank.resize( 1.5 * _batch_size_initial );

  //populate initial source bank
  for ( auto i = 0 ; i < _batch_size_initial ; ++i ) {
    const auto starting_cell = geom::find_cell( pos, cells );
    _source_bank[i] = mc::fissionSite( pos, 1.5, starting_cell );
  }

  if ( _dos_estimators    ) _dos_estimators->resize_source_bank(    1.5 * _batch_size_initial );
  if ( _dos_ur_estimators ) _dos_ur_estimators->resize_source_bank( 1.5 * _batch_size_initial );
}

mc::particle mc::kEigenvalueCalc::source_particle( const uint32_t id, util::RNG::ptr rng ) const {
  const auto mu  = 2.0*rng->sample() - 1.0;
  const auto azi = constants::two_pi * rng->sample();
  const auto sin_theta = std::sqrt( 1.0 - mu*mu );

  particle p;
  p.pos  = _source_bank[id].pos;
  p.dir  = util::point( mu, sin_theta*std::cos(azi), sin_theta*std::sin(azi) );
  p.erg  = _source_bank[id].erg;
  p.wgt  = _source_weight;
  p.cell =  _source_bank[id].cell_ptr;

  if ( _dos_estimators    ) _dos_estimators->start_history(    _source_bank[id].dos_id );
  if ( _dos_ur_estimators ) _dos_ur_estimators->start_history( _source_bank[id].dos_id );

  return p;
}

void mc::kEigenvalueCalc::score( const estimator_type_t type, const double contribution, const uint32_t geom_id,
                                 const mc::particle& p, const xs::xsBuffer& xs_buffer, util::RNG::ptr rng ) {

  if ( type == estimator_type_t::TRACK_LENGTH_FLUX ) {
    //score internal keff track-length estimator
    _keff_est += xs_buffer.nufission_macro_xs * contribution;

    //note: values cached in score_indirect_stream and used in score_direct_responses
    //      therefore cannot switch the order of these calls
    if ( _dos_estimators ) {
      _dos_estimators->score_indirect_stream( contribution, p.erg, p.cell, xs_buffer );
      if ( is_scoring_batch() ) {
        _dos_estimators->score_direct_response( contribution, p.erg, p.cell, xs_buffer );
      }
    }
    if ( _dos_ur_estimators ) {
      _dos_ur_estimators->score_indirect_stream( contribution, p.erg, p.cell, xs_buffer );
      if ( is_scoring_batch() ) {
        _dos_ur_estimators->score_direct_response( contribution, p.erg, p.cell, xs_buffer );
      }
    }
  }
  else if ( type == estimator_type_t::COLLISION_FLUX ) {
    //bank fission neutrons at each collision
    const auto num_fission_neutrons = std::floor( p.wgt * xs_buffer.nufission_macro_xs / ( _keff_norm * xs_buffer.total_macro_xs ) + rng->sample() );

    //compute differential operator sampling fission source contributions if present
    uint32_t dos_id = 0;
    if ( _dos_estimators ) {
      if ( num_fission_neutrons > 0 ) {
        dos_id = _dos_estimators->insert_fission_source( p.wgt, p.erg, p.cell, num_fission_neutrons, xs_buffer );
      }
      _dos_estimators->score_indirect_collision( p.wgt, p.erg, p.cell, xs_buffer );
    }
    if ( _dos_ur_estimators ) {
      if ( num_fission_neutrons > 0 ) {
        dos_id = _dos_ur_estimators->insert_fission_source( p.wgt, p.erg, p.cell, num_fission_neutrons, xs_buffer );
      }
      _dos_ur_estimators->score_indirect_collision( p.wgt, p.erg, p.cell, xs_buffer );
    }

    //bank the new source neutrons
    const auto mat = p.cell->material();
    for ( auto k = 0 ; k < num_fission_neutrons ; ++k ) {
      //sample fission emitting isotope for each neutron
      const auto iso = mat->select_isotope( xs_buffer.nufission_macro_xs, xs_buffer.nufission_xs, rng );

      //sample outoing neutron energy
      auto [ Eout, dirout ] = iso->sample_fission( p.erg, p.dir, xs_buffer, rng );

      //bank fission neutron
      _fission_bank[_nbanked] = mc::fissionSite( p.pos, Eout, p.cell );
      _fission_bank[_nbanked].dos_id = dos_id;
      ++_nbanked;
    }
  }
  if ( is_scoring_batch() ) {
    score_estimators_of_type( type, _estimators, contribution, geom_id, p.erg, p.cell, xs_buffer );
  }
}

void mc::kEigenvalueCalc::end_history() {
  if (  is_scoring_batch() ) {
    //sensitivity coefficients must be scored before estimators are cleared
    if ( _dos_estimators    ) _dos_estimators->end_history();
    if ( _dos_ur_estimators ) _dos_ur_estimators->end_history();

    end_history_estimators( _estimators );
  }
}

void mc::kEigenvalueCalc::end_batch() {
  if ( _dos_estimators    ) _dos_estimators->end_batch( _nbanked );
  if ( _dos_ur_estimators ) _dos_ur_estimators->end_batch( _nbanked );

  _keff_est /= _batch_size_initial;
  _source_weight = static_cast<double>( _batch_size_initial ) / _nbanked;
  _keff_norm     = _keff_est;
  _batch_size    = _nbanked;

  std::cout << std::fixed << std::setprecision(8);
  if ( _batch >= _num_inactive ) {
    _mean_keff += _keff_est;
    _stdv_keff += _keff_est * _keff_est;

    const auto nactive = _batch - _num_inactive + 1;
    const auto print_mean_keff = _mean_keff / nactive;
    const auto print_stdv_keff =
      nactive > 1 ?
      std::sqrt( ( _stdv_keff/nactive - print_mean_keff * print_mean_keff )/(nactive-1) ) : 0.0;
    std::cout << std::setw(5) << _batch+1 << ' ' << _keff_est << ' ' << print_mean_keff << ' ' << print_stdv_keff << '\n';

//    if ( _dos_estimators    ) _dos_estimators->write_sensitivity( static_cast<double>( _batch_size_initial ) * nactive );
//    if ( _dos_ur_estimators ) _dos_ur_estimators->write_sensitivity( static_cast<double>( _batch_size_initial ) * nactive );
  }
  else {
    std::cout << std::setw(5) << _batch+1 << ' ' << _keff_est << ' ' << _nbanked << '\n';
  }

  std::swap( _source_bank, _fission_bank );
  _keff_est = 0.0;
  _nbanked  = 0;
  _batch++;

  if ( is_done() ) {
    write_estimators( _estimators, normalization() );
    if ( _dos_estimators    ) _dos_estimators->write_sensitivity(    normalization() );
    if ( _dos_ur_estimators ) _dos_ur_estimators->write_sensitivity( normalization() );
  }
}
