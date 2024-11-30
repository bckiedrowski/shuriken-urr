#include <cmath>
#include <iomanip>
#include <iostream>

#include "constants.hpp"

#include "fixedSourceCalc.hpp"

mc::particle mc::fixedSourceCalc::source_particle( const uint32_t id, util::RNG::ptr rng ) const {
  if ( _dos_estimators    ) _dos_estimators->start_history( id );
  if ( _dos_ur_estimators ) _dos_ur_estimators->start_history( id );
  return _src_dist( rng );
}

void mc::fixedSourceCalc::score( const estimator_type_t type, const double contribution, const uint32_t geom_id,
                                 const mc::particle& p, const xs::xsBuffer& xs_buffer, util::RNG::ptr rng ) {
  if ( type == estimator_type_t::TRACK_LENGTH_FLUX ) {
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
    if ( _dos_estimators    ) _dos_estimators->score_indirect_collision(    p.wgt, p.erg, p.cell, xs_buffer );
    if ( _dos_ur_estimators ) _dos_ur_estimators->score_indirect_collision( p.wgt, p.erg, p.cell, xs_buffer );
  }
  if (  is_scoring_batch() ) {
    score_estimators_of_type( type, _estimators, contribution, geom_id, p.erg, p.cell, xs_buffer );
  }
}

void mc::fixedSourceCalc::end_history() {
  if (  is_scoring_batch() ) {
    //sensitivity coefficients must be scored before estimators are cleared
    if ( _dos_estimators    ) _dos_estimators->end_history();
    if ( _dos_ur_estimators ) _dos_ur_estimators->end_history();

    end_history_estimators( _estimators );
  }
}

void mc::fixedSourceCalc::end_batch() {
  _batch++;
  std::cout << '\n' << " finished bactch " << _batch << " with " << _batch * _batch_size << " histories completed" << '\n';

  write_estimators( _estimators, normalization() );
  if ( _dos_estimators    ) _dos_estimators->write_sensitivity(    normalization() );
  if ( _dos_ur_estimators ) _dos_ur_estimators->write_sensitivity( normalization() );
}
