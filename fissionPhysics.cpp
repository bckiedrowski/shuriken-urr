#include "constants.hpp"

#include "fissionPhysics.hpp"

xs::fissionPhysics::fissionPhysics( const aceData& acefile, const uint32_t mt_offset ) :
  scatterPhysics( acefile, mt_offset ) {
  _num_delayed_groups = acefile.nxs(8);

  //note: if this is zero, then only prompt fission data is available and nothing else
  //      needs to be done here since prompt data is constructed in the base constructor
  if ( _num_delayed_groups > 0 ) {
    const auto loc_nu     = acefile.jxs(2);  //location of fission yields
    const auto loc_dnu    = acefile.jxs(24); //location of delayed yield data
    const auto loc_bdd    = acefile.jxs(25); //location of basic delayed data
    const auto loc_dnedl  = acefile.jxs(26); //location of delayed neutron spectra locators
    const auto loc_dned   = acefile.jxs(27); //location of delayed neutron spectra data

    if ( acefile.xss( loc_nu ) > 0 && loc_dnu > 0 ) {
      //prompt (or total) and delayed yield data are given
      //TODO: implement this case
      assert( false );
    }
    else if ( acefile.xss( loc_nu ) < 0 ) {
      //both prompt and total yield data are given
      const auto nu_offset = std::abs( static_cast<int32_t>( acefile.xss( loc_nu ) ) );

      const auto loc_prompt_nu = loc_nu + 1;
      const auto loc_total_nu  = loc_nu + nu_offset + 1;

      const auto prompt_yield_type = static_cast<int32_t>( acefile.xss(loc_prompt_nu) );
      const auto total_yield_type  = static_cast<int32_t>( acefile.xss(loc_total_nu)  );

      //create prompt yield data based on type
      if ( prompt_yield_type == 1 ) {
        _prompt_yield = std::make_shared<polynomialYield>( acefile, loc_prompt_nu+1 );
      }
      else if ( prompt_yield_type == 2 ) {
        _prompt_yield = std::make_shared<tabularYield>( acefile, loc_prompt_nu+1 );
      }
      else {
        assert("unknown fission yield type flag.");
      }

      //create total yield data based on type
      if ( total_yield_type == 1 ) {
        _total_yield = std::make_shared<polynomialYield>( acefile, loc_total_nu+1 );
      }
      else if ( prompt_yield_type == 2 ) {
        _total_yield = std::make_shared<tabularYield>( acefile, loc_total_nu+1 );
      }
      else {
        assert("unknown fission yield type flag.");
      }
    }
    else {
      //no yield data present, file likely has an error
      assert(false);
    }

    //construct delayed neutron group decay constant and probability vectors
    {
      auto loc = loc_bdd;
      _decay_constants.resize( _num_delayed_groups );
      _delayed_group_prob.resize( _num_delayed_groups );
      for ( auto i = 0 ; i < _num_delayed_groups ; ++i ) {
        _decay_constants[i]    = acefile.xss( loc ) * constants::shakes_per_second;
        _delayed_group_prob[i] = std::make_shared< interpGrid >( acefile, loc+1 );
        //advance locator to the next group
        const auto nr = static_cast<uint32_t>( acefile.xss(loc+1) );
        const auto ne = static_cast<uint32_t>( acefile.xss(loc+2*nr+2) );
        loc += 2*(nr+ne) + 3;
      }
    }

    //construct inelastic scattering laws for delayed neutron emission spectra
    {
      _group_law_offsets.resize( _num_delayed_groups + 2 );
      _group_law_offsets[0] = 0;
      for ( auto i = 1 ; i <= _num_delayed_groups ; ++i ) {
        _group_law_offsets[i] = _laws.size();
        create_laws( acefile, loc_dnedl, loc_dned, i );
      }
      _group_law_offsets[ _num_delayed_groups + 1 ] = _laws.size();
    }

  }
}

std::pair<double,util::point> xs::fissionPhysics::sample(
  const double erg, const util::point& dir, util::RNG::ptr rng ) const {
  //determine the delayed group
  auto igroup = 0;
  if ( _num_delayed_groups > 0 ) {
    const auto r1 = rng->sample();
    //check prompt fission first
    const auto p_prompt = _prompt_yield->yield( erg ) / _total_yield->yield( erg );
    if ( r1 < p_prompt ) {
      //assign prompt neutron
      igroup = 0;
    }
    else {
      auto r2 = rng->sample();
      auto ig = 0;
      //note if second to last group not selected, takes last group automatically
      //this avoids issues with roundoff
      for ( ; ig < _num_delayed_groups-1 ; ++ig ) {
        r2 -= _delayed_group_prob[ig]->operator()( erg );
        if ( r2 < 0 ) break;
      }
      igroup = ig+1;
    }
  }
  else {
    //assign prompt neutron
    igroup = 0;
  }

  //sample the scattering law (if there are multiple) based on prescribed probability
  auto ilaw = _group_law_offsets[igroup];
  if ( _group_law_offsets[igroup+1] - ilaw > 1 ) {
    auto r = rng->sample();
    for ( ; ilaw < _group_law_offsets[igroup+1]-1 ; ++ilaw ) {
      //reduce random number by probability and take once falls below zero
      //note loop index is such that the final law is taken otherwise (avoids roundoff issues)
      r -= _laws[ilaw].prob->operator()( erg );
      if ( r <= 0.0 ) break;
    }
  }
  return _laws[ilaw].data->sample( erg, dir, rng );
}
