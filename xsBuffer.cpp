#include "xsBuffer.hpp"


xs::xsBuffer::xsBuffer( const uint32_t niso ) {
  erg.resize( niso, 0.0 );
  total_xs.resize( niso, 0.0 );
  elastic_xs.resize( niso, 0.0 );
  fission_xs.resize( niso, 0.0 );
  nufission_xs.resize( niso, 0.0 );
  ngamma_xs.resize( niso, 0.0 );
  urr_rn.resize( niso, 0.0 );

  urr_data.resize( niso );
}

void xs::xsBuffer::flush() {
  for ( auto i = 0 ; i < erg.size() ; ++i ) {
    erg[i]    = 0.0;
  }
}

double xs::xsBuffer::operator()( const xs::endf_mt_t mt, const uint32_t id ) const {
  switch( mt ) {
  case xs::endf_mt_t::TOTAL:
    return total_xs[id];
  case xs::endf_mt_t::ELASTIC:
    return elastic_xs[id];
  case xs::endf_mt_t::FISSION:
    return fission_xs[id];
  case xs::endf_mt_t::Z_GAMMA:
    return ngamma_xs[id];
  case xs::endf_mt_t::NUFISSION:
    return nufission_xs[id];
  default:
    //error to attempt an access of nonexistent element
    assert( false );
    return 0.0;
  }
}

std::pair< double, double > xs::xsBuffer::urr_table_xs( const xs::endf_mt_t mt, const uint32_t id ) const {
  switch( mt ) {
  case xs::endf_mt_t::TOTAL:
    return std::make_pair( urr_data[id].left_elastic_xs  + urr_data[id].left_fission_xs  + urr_data[id].left_ngamma_xs,
                           urr_data[id].right_elastic_xs + urr_data[id].right_fission_xs + urr_data[id].right_ngamma_xs );
  case xs::endf_mt_t::ELASTIC:
    return std::make_pair( urr_data[id].left_elastic_xs, urr_data[id].right_elastic_xs );
  case xs::endf_mt_t::FISSION:
    return std::make_pair( urr_data[id].left_fission_xs, urr_data[id].right_fission_xs );
  case xs::endf_mt_t::Z_GAMMA:
    return std::make_pair( urr_data[id].left_ngamma_xs,  urr_data[id].right_ngamma_xs  );
  case xs::endf_mt_t::NUFISSION:
    return std::make_pair( urr_data[id].left_nufission_xs, urr_data[id].right_nufission_xs );
  default:
    //error to attempt an access of nonexistent element
    assert( false );
    return std::make_pair( 0.0, 0.0 );
  }
}
