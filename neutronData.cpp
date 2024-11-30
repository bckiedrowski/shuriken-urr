#include <cmath>
#include <iostream>

#include "binarySearch.hpp"
#include "neutronData.hpp"

bool xs::neutronData::read_ace( const aceData& acefile ) {

  _zaid = acefile.zaid();
  _temp = acefile.temp();

  _nerg = acefile.nxs(3);
  _nrxn = 1 + acefile.nxs(4);  //nxs(4) excludes elastic

  _egrid.resize( _nerg );
  _mt.resize( _nrxn );

  //see page F-4 of mcnp5 manual vol. III for detailed definitions of these blocks
  const auto loc_esg    = acefile.jxs(1);
  const auto loc_nu     = acefile.jxs(2);
  const auto loc_mtr    = acefile.jxs(3);
  const auto loc_lqr    = acefile.jxs(4);
  const auto loc_tyr    = acefile.jxs(5);
  const auto loc_lsig   = acefile.jxs(6);  //location of individual xs table locators
  const auto loc_sig    = acefile.jxs(7);  //location of xs (SIG) block in xss
  const auto loc_land   = acefile.jxs(8);  //location of angular distribution locators
  const auto loc_and    = acefile.jxs(9);  //location of angular distribution block
  const auto loc_dlw    = acefile.jxs(11);

  const auto nang_rxn = 1 + acefile.nxs(5);  //number of reactions with angular distributions

  //elastic scattering is special in ace format
  //copy into map along with all other reactions for simplicity
  const auto loc_elastic = loc_esg + _nerg*3;
  rxn_t elastic;
  elastic.xs.resize( _nerg );
  elastic.ierg   = 0;     // follow C index convention, not Fortran
  elastic.nerg   = _nerg;
  elastic.qvalue = 0.0;
  elastic.frame  = ref_frame_t::COM_FRAME;
  elastic.yield  = std::make_shared<fixedYield>(1);

  //populate energy grid and elastic cross sections
  for ( auto i = 0 ; i < _nerg ; ++i ) {
    _egrid[i]     = acefile.xss( loc_esg   + i );
    elastic.xs[i] = acefile.xss( loc_elastic + i );
  }
  //scattering physics
  elastic.scatter_physics = std::make_shared< scatterPhysics >( acefile, 0 );

  //insert elastic (MT=2) into reaction map
  _rxn.insert( std::pair< endf_mt_t, rxn_t >( endf_mt_t::ELASTIC, elastic ) );
  _mt[0] = endf_mt_t::ELASTIC;

  //read the other reactions (-1 offset because elastic is treated special)
  for ( auto i = 0 ; i < _nrxn-1 ; ++i ) {
    rxn_t this_rxn;

    _mt[i+1] = static_cast<endf_mt_t>( acefile.xss( loc_mtr + i ) );

    const auto ty  = static_cast<int32_t>( acefile.xss( loc_tyr + i ) );
    if ( std::abs( ty ) <= 4 ) {
      //see pg. F-16 in mcnp5 manual for explanation of constants
      //fixed yield
      this_rxn.yield = std::make_shared<fixedYield>( std::abs(ty) );
    }
    else if ( std::abs( ty ) == 19 ) {
      //store energy-dependent TOTAL yield in yield structure for now
      if ( loc_nu > 0 ) {
        if ( acefile.xss(loc_nu) > 0 ) {
          //only prompt or total data given, treat as total for now
          //note: offset locator by 1 to be start of the table beyond flags
          const auto yield_type_flag = static_cast<int32_t>( acefile.xss(loc_nu) );
          if ( yield_type_flag == 1 ) {
            this_rxn.yield = std::make_shared<polynomialYield>( acefile, loc_nu+1 );
          }
          else if ( yield_type_flag == 2 ) {
            this_rxn.yield = std::make_shared<tabularYield>( acefile, loc_nu+1 );
          }
          else {
            assert("unknown fission yield type flag.");
          }
        }
        else {
          //both prompt and total data given, store only total here
          const auto loc_total_nu = loc_nu + std::abs( static_cast<int32_t>(acefile.xss(loc_nu)) ) + 1;
          const auto yield_type_flag = static_cast<int32_t>( acefile.xss(loc_total_nu) );
          if ( yield_type_flag == 1 ) {
            this_rxn.yield = std::make_shared<polynomialYield>( acefile, loc_total_nu+1 );
          }
          else if ( yield_type_flag == 2 ) {
            this_rxn.yield = std::make_shared<tabularYield>( acefile, loc_total_nu+1 );
          }
          else {
            assert("unknown fission yield type flag.");
          }
        }
     }
     else {
       assert( false );
     }

    }
    else if ( std::abs( ty ) > 100 ) {
      //energy-dependent tabular neutron yield
      const auto loc_yield = loc_dlw + std::abs(ty) - 101;
      this_rxn.yield = std::make_shared<tabularYield>( acefile, loc_yield );
    }
    else {
      //undefined yield number
      assert( false );
    }
    this_rxn.frame  = ty < 0 ? ref_frame_t::COM_FRAME : ref_frame_t::LAB_FRAME;
    this_rxn.qvalue = acefile.xss( loc_lqr + i );

    const auto loca = static_cast<int32_t>( acefile.xss( loc_lsig + i ) );
    this_rxn.ierg   = static_cast<int32_t>( acefile.xss( loc_sig + loca - 1 ) ) - 1; //follow C array indexing
    this_rxn.nerg   = static_cast<int32_t>( acefile.xss( loc_sig + loca     ) );

    this_rxn.xs.resize( this_rxn.nerg );
    for ( auto k = 0 ; k < this_rxn.nerg ; ++k ) {
      this_rxn.xs[k] = acefile.xss( loc_sig + loca + 1 + k );
    }

    //energy-angle scattering distributions
    //reactions ordering in array is order such that the ones with scattering distributions are first
    if ( i < nang_rxn-1 ) {
      //fission is handled specially because of delayed neutrons
      const auto mt = static_cast<endf_mt_t>(_mt[i+1]);
      if ( mt == endf_mt_t::FISSION || mt == endf_mt_t::Z_F || mt == endf_mt_t::Z_NF
        || mt == endf_mt_t::Z_2NF   || mt == endf_mt_t::Z_3NF ) {
        this_rxn.scatter_physics = std::make_shared< fissionPhysics >( acefile, i+1 );
      }
      else {
        this_rxn.scatter_physics = std::make_shared< scatterPhysics >( acefile, i+1 );
      }
    }
    _rxn.insert( std::pair< endf_mt_t, rxn_t >( _mt[i+1], this_rxn ) );
  }
  //insert dummy to terminate xs search
  _rxn.emplace( std::pair< endf_mt_t, rxn_t >( endf_mt_t::TERMINATE, rxn_t() ) );

  //insert heating data into reaction array
  rxn_t heating;
  heating.xs.resize( _nerg, 0.0 );
  heating.ierg   = 0;     // follow C index convention, not Fortran
  heating.nerg   = _nerg;
  heating.qvalue = 0.0;
  heating.frame  = ref_frame_t::UNKNOWN;

  const auto loc_heating = loc_esg + _nerg*4;
  for ( auto i = 0 ; i < _nerg ; ++i ) {
    heating.xs[i] = acefile.xss( loc_heating + i );
  }
  _rxn.insert( std::pair< endf_mt_t, rxn_t >( endf_mt_t::HEATING, heating ) );

  //compute total cross section
  rxn_t total;
  total.xs.resize( _nerg, 0.0 );
  total.ierg   = 0;     // follow C index convention, not Fortran
  total.nerg   = _nerg;
  total.qvalue = 0.0;
  total.frame  = ref_frame_t::UNKNOWN;

  bool mt18_fission_found = false;
  for ( auto it = _rxn.begin() ; it != _rxn.end() ; ++it ) {
    const auto mt = it->first;

    //exit loop if termination point hit
    if ( mt == endf_mt_t::TERMINATE ) break;

    //skip inelastic if found, accumulate partials
    if ( mt == endf_mt_t::INELASTIC ) continue;

    //special case for fission, if MT=18 given, do not read partials MT=19,20,21,38
    if ( mt == endf_mt_t::FISSION ) mt18_fission_found = true;
    if ( mt18_fission_found
         && ( mt == endf_mt_t::Z_F   || mt == endf_mt_t::Z_NF
           || mt == endf_mt_t::Z_2NF || mt == endf_mt_t::Z_3NF ) ) continue;

    //record MT=18 not given and record it as given as partials or chances
    if ( mt == endf_mt_t::Z_F ) _fission_chances_given = true;

    // populate total cross section on energy grid
    const auto& rxn   = it->second;
    const auto  start = rxn.ierg;
    for ( auto k = 0 ; k < rxn.nerg ; ++k ) {
      total.xs[start+k] += rxn.xs[k];
    }
  }
  _rxn.insert( std::pair< endf_mt_t, rxn_t >( endf_mt_t::TOTAL, total ) );

  //set fissionable flag
  _fissionable = mt18_fission_found || _fission_chances_given;

  //unresolved resonance data
  if ( acefile.jxs(23) > 0 ) {
    _urr = std::make_shared< urrData >( acefile );
//    _urr->unionize();
    _urr_exist = true;
  }

  //set maximum energy in case input energy runs off the table
  _emax = _egrid[ _egrid.size() - 1 ] - 1.0e-11;

  //CANNOT CHANGE _rxn AFTER THIS POINT OR ITERATORS WILL BE INVALIDATED
  //set special iterators to accelerate calculations and handle unresolved resonances
  //note that these reactions are usually the most common and are checked first
  //when sampling reaction type
  _rxn_total_iter     = _rxn.find( endf_mt_t::TOTAL     );
  _rxn_elastic_iter   = _rxn.find( endf_mt_t::ELASTIC   );
  _rxn_ngamma_iter    = _rxn.find( endf_mt_t::Z_GAMMA   );
  _rxn_terminate_iter = _rxn.find( endf_mt_t::TERMINATE );
  if ( _fissionable ) {
    //need to handle two possible cases for handling fission
    if ( _fission_chances_given ) {
      _rxn_fission_iter = _rxn.find( endf_mt_t::Z_F );
      _rxn_fission_chance_iter[0] =  _rxn.find( endf_mt_t::Z_NF  );
      _rxn_fission_chance_iter[1] =  _rxn.find( endf_mt_t::Z_2NF );
      _rxn_fission_chance_iter[2] =  _rxn.find( endf_mt_t::Z_3NF );
    }
    else {
     _rxn_fission_iter = _rxn.find( endf_mt_t::FISSION );
    }
  }

  //build vector of other reactions (usually these are minor contributions)
  for ( auto it = _rxn_elastic_iter ; it != _rxn_terminate_iter ; ++it ) {
    const auto mt = it->first;

    //do not take any of the following
    if ( mt == endf_mt_t::ELASTIC   ) continue;
    if ( mt == endf_mt_t::INELASTIC ) continue;
    if ( mt == endf_mt_t::FISSION   ) continue;
    if ( mt == endf_mt_t::Z_F       ) continue;
    if ( mt == endf_mt_t::Z_NF      ) continue;
    if ( mt == endf_mt_t::Z_2NF     ) continue;
    if ( mt == endf_mt_t::Z_3NF     ) continue;
    if ( mt == endf_mt_t::Z_GAMMA   ) continue;

    _other_rxn_iter.push_back( it );
  }

  //create hash table to limit time spent in binary search
  _hash_min   = erg_hash_fxn( _egrid[0] );
  _hash_max   = erg_hash_fxn( _egrid[_egrid.size()-1] );
  _hash_delta = nhash/( _hash_max - _hash_min );
  auto ih = 1;
  _hash_lookup[0] = 0;
  auto hash_next = 1.0/_hash_delta;

  for ( auto ie = 1 ; ie < _egrid.size()-1 ; ++ie ) {
    const auto h = erg_hash_fxn( _egrid[ie] ) - _hash_min;
    while ( h > hash_next ) {
      _hash_lookup[ih] = ie-1;
      hash_next += 1.0/_hash_delta;
      ++ih;
    }
  }
  // fill in remainder
  for ( ; ih <= nhash ; ++ih ) {
    _hash_lookup[ih] = _egrid.size() - 2;
  }

  return true; //success
}


double xs::neutronData::erg_hash_fxn( const double x ) const {
  int b;
  const auto a = std::frexp( x, &b );
  return (b-1) + 2.0 * ( a - 0.5 );
};

std::tuple<uint32_t,uint32_t,uint32_t> xs::neutronData::split_zaid() const {
  //Z (charge number), A (mass number), S (state number)
  return std::make_tuple( ( _zaid / 1000 ) % 1000 , _zaid % 1000, _zaid / 1000000 );
}

//utility function to do find a cross section
//this routine is more expensive than fast_xs_lookup and should only be used when
//a single isolated cross section is needed or when the iterator to the rxn is
//not easily available
double xs::neutronData::find_xs( const endf_mt_t mt, const double erg ) const {
  const auto& rxn   = _rxn.at(mt);
  const auto  start = rxn.ierg;  //start of the portion of energy grid valid for rxn
  const auto  end   = start + rxn.nerg - 1;  //end of portion of energy grid

  const auto& xs_data = rxn.xs;

  //check out of bounds, return value endpoint
  if ( erg <  _egrid[start] ) { return xs_data[0]; }
  if ( erg >= _egrid[end]   ) { return xs_data[rxn.nerg-1]; }

  //binary search on global energy grid fine grid of cross section
  const auto ih  = std::floor( _hash_delta*( erg_hash_fxn(erg) - _hash_min ) );
  const auto ie  = util::binary_search( erg, _egrid, _hash_lookup[ih], 1 + _hash_lookup[ih+1] );
  const auto rtp = ( erg - _egrid[ie] )/( _egrid[ie+1] - _egrid[ie] );

  const auto k = ie - start;
  return xs_data[k] + rtp*( xs_data[k+1] - xs_data[k] );
}

//grab cross section on provided index of global energy grid k and interpolation fraction r
double xs::neutronData::fast_xs_lookup( const rxn_t& rxn, const uint32_t k, const double r ) const {
  const auto  start = rxn.ierg;  //start of the portion of energy grid valid for rxn
  const auto  end   = start + rxn.nerg - 1;  //end of portion of energy grid

  const auto& xs_data = rxn.xs;

  //return lower and upper grid points if off edge of table
  if ( k <  start ) { return xs_data[0]; }
  if ( k >= end   ) { return xs_data[rxn.nerg-1]; }

  //within table, do lin-lin interpolation
  const auto i = k - start;
  return xs_data[i] + r*( xs_data[i+1] - xs_data[i] );
}

void xs::neutronData::fill_buffer( const double erg, xs::xsBuffer& xs_buffer, util::RNG::ptr rng ) const {
  const auto k = _id;
  xs_buffer.erg[k] = erg;

  if ( ! within_urr( erg ) ) {
    const auto E   = std::fmin( erg, _emax );
    const auto ih  = std::floor( _hash_delta*( erg_hash_fxn(E) - _hash_min ) );
    const auto ie  = util::binary_search( E, _egrid, _hash_lookup[ih], 1 + _hash_lookup[ih+1] );
    const auto rtp = ( E - _egrid[ie] )/( _egrid[ie+1] - _egrid[ie] );

    xs_buffer.total_xs[k]   = fast_xs_lookup( _rxn_total_iter->second,   ie, rtp );
    xs_buffer.elastic_xs[k] = fast_xs_lookup( _rxn_elastic_iter->second, ie, rtp );
    xs_buffer.ngamma_xs[k]  = fast_xs_lookup( _rxn_ngamma_iter->second,  ie, rtp );

    if ( _fissionable ) {
      const auto nubar = _rxn_fission_iter->second.yield->yield( E );

      //TODO: TESTING CODE -- REMOVE THIS AND UNCOMMENT CODE BELOW!
      //const auto factor = ( _zaid == 92235 && erg >= 0.1 && erg < 0.5 ) ? 1.01 : 1.0;
      //const auto fxs_base  = fast_xs_lookup( _rxn_fission_iter->second, ie, rtp );
      //const auto fxs_pert  = factor * fxs_base;
      //const auto fxs_delta = fxs_pert - fxs_base;

      //xs_buffer.fission_xs[k]    = fxs_pert;
      //xs_buffer.nufission_xs[k]  = nubar * xs_buffer.fission_xs[k];
      //xs_buffer.total_xs[k]     += fxs_delta;

      xs_buffer.fission_xs[k]   = fast_xs_lookup( _rxn_fission_iter->second, ie, rtp );
      xs_buffer.nufission_xs[k] = nubar * xs_buffer.fission_xs[k];
      if ( _fission_chances_given ) {
        //if fission given as chances/partials, then first chance was already handled,
        //but need to consider other chances with thresholds of increasing energy
        //and add them as appropriate
        //* note that usually the yield (nu) for each fission chance is the same in the
        //  data, so this could be accelerated. it is theoretically possible this could
        //  change in the future, so leaving the code like this in a slightly less
        //  efficient form
        for ( auto it : _rxn_fission_chance_iter ) {
          const auto fxs = fast_xs_lookup( it->second, ie, rtp );
          if ( fxs == 0.0 ) break;
          const auto nu  = it->second.yield->yield( E );
          xs_buffer.fission_xs[k]   += fxs;
          xs_buffer.nufission_xs[k] += nu*fxs;
        }
      }
    }
    else {
      xs_buffer.fission_xs[k]   = 0.0;
      xs_buffer.nufission_xs[k] = 0.0;
    }
  }
  else {
    const auto rn = rng->sample();
    const auto urr_xs = sample_urr( erg, rn );

    xs_buffer.urr_rn[_id]     = rn;
    xs_buffer.total_xs[_id]   = urr_xs.total;
    xs_buffer.elastic_xs[_id] = urr_xs.elastic;
    xs_buffer.fission_xs[_id] = urr_xs.fission;
    xs_buffer.ngamma_xs[_id]  = urr_xs.ngamma;

    //urr data structure for dos sensitivities
    xs_buffer.urr_data[_id].energy_index    = urr_xs.energy_index;
    xs_buffer.urr_data[_id].interp_factor   = urr_xs.interp_factor;

    xs_buffer.urr_data[_id].left_xs_index   = urr_xs.left_xs_index;
    xs_buffer.urr_data[_id].left_elastic_xs = urr_xs.left_elastic_xs;
    xs_buffer.urr_data[_id].left_fission_xs = urr_xs.left_fission_xs;
    xs_buffer.urr_data[_id].left_ngamma_xs  = urr_xs.left_ngamma_xs;

    xs_buffer.urr_data[_id].right_xs_index   = urr_xs.right_xs_index;
    xs_buffer.urr_data[_id].right_elastic_xs = urr_xs.right_elastic_xs;
    xs_buffer.urr_data[_id].right_fission_xs = urr_xs.right_fission_xs;
    xs_buffer.urr_data[_id].right_ngamma_xs  = urr_xs.right_ngamma_xs;

    if ( _fissionable ) {
      //note higher chances of fission do not occur in unresolved region
      //because they would normally be well above the resonance region
      const auto nubar = _rxn_fission_iter->second.yield->yield( erg );
      xs_buffer.nufission_xs[k] = nubar * xs_buffer.fission_xs[k];

      xs_buffer.urr_data[_id].left_nufission_xs  = nubar * urr_xs.left_fission_xs;
      xs_buffer.urr_data[_id].right_nufission_xs = nubar * urr_xs.right_fission_xs;
    }
  }
};

xs::endf_mt_t xs::neutronData::select_rxn( const double erg, const xs::xsBuffer& xs_buffer, util::RNG::ptr rng ) const {
  //assumes xs_buffer has been populated by fill_buffer routine
  assert( xs_buffer.erg[_id] == erg );

  const auto k = _id;
  auto t = xs_buffer.total_xs[k] * rng->sample();

  //explicitly go through major reactions (which are also on unresolved table)
  t -= xs_buffer.elastic_xs[k];
  if ( t < 0 ) return endf_mt_t::ELASTIC;

  t -= xs_buffer.ngamma_xs[k];
  if ( t < 0 ) return endf_mt_t::Z_GAMMA;

  t -= xs_buffer.fission_xs[k];
  if ( t < 0 ) return endf_mt_t::FISSION;

  //if not found, loop over the other reactions
  const auto E   = std::fmin( erg, _emax );
  const auto ih  = std::floor( _hash_delta*( erg_hash_fxn(E) - _hash_min ) );
  const auto ie  = util::binary_search( E, _egrid, _hash_lookup[ih], 1 + _hash_lookup[ih+1] );
  const auto rtp = ( E - _egrid[ie] )/( _egrid[ie+1] - _egrid[ie] );
  for ( auto it : _other_rxn_iter ) {
    t -= fast_xs_lookup( it->second, ie, rtp );
    if ( t < 0 ) return it->first;
  }
  //if this point is reached, the total and partials do not balance
  //print out results for debugging purposes and return an unknown type
  auto total = xs_buffer.elastic_xs[k] + xs_buffer.ngamma_xs[k] + xs_buffer.fission_xs[k];
  for ( auto it : _other_rxn_iter ) {
    std::cout << static_cast<int>( it->first ) << ' ' << fast_xs_lookup( it->second, ie, rtp ) << '\n';
    total += fast_xs_lookup( it->second, ie, rtp );
  }
  std::cout << total << ' ' << xs_buffer.total_xs[k] << '\n';

  assert( false );
  return endf_mt_t::NONE;
}

double xs::neutronData::yield( const endf_mt_t mt, const double erg ) const {
  return _rxn.at(mt).yield->yield( erg );
}

bool xs::neutronData::within_urr( const double erg ) const {
  if ( _urr_exist ) {
    return erg >= _urr->lower_bound() && erg < _urr->upper_bound();
  }
  //no unresolved probability tables for this isotope
  return false;
}

xs::urrXS xs::neutronData::sample_urr( const double erg, const double random_number ) const {
  return _urr->xs( erg, random_number );
}

std::pair<double,util::point> xs::neutronData::sample_secondary(
  const endf_mt_t mt, const double erg, const util::point& dir, util::RNG::ptr rng ) const {
  //return energy and direction vector pair
  return _rxn.at(mt).scatter_physics->sample( erg, dir, rng );
}

std::pair<double,util::point> xs::neutronData::sample_fission(
  const double erg, const util::point& dir, const xs::xsBuffer& xs_buffer, util::RNG::ptr rng ) const {
  assert( _fissionable );

  const auto E = std::fmin( erg, _emax );
  if ( _fission_chances_given && (! within_urr( erg )) ) {
    //sample chance of fission based on nu fission cross section ratios
    //and select appropriate scattering law
    const auto ih = std::floor( _hash_delta*( erg_hash_fxn(E) - _hash_min ) );
    const auto ie = util::binary_search( E, _egrid, _hash_lookup[ih], 1 + _hash_lookup[ih+1] );
    const auto rtp = ( E - _egrid[ie] )/( _egrid[ie+1] - _egrid[ie] );

    auto t = xs_buffer.nufission_xs[_id] * rng->sample();
    t -= _rxn_fission_iter->second.yield->yield( E ) * fast_xs_lookup( _rxn_fission_iter->second, ie, rtp );
    if ( t < 0.0 ) return _rxn_fission_iter->second.scatter_physics->sample( E, dir, rng );
    for ( auto it : _rxn_fission_chance_iter ) {
      t -= it->second.yield->yield( E ) * fast_xs_lookup( it->second, ie, rtp );
      if ( t < 0.0 ) return it->second.scatter_physics->sample( E, dir, rng );
    }
    assert( false );
    return _rxn_fission_iter->second.scatter_physics->sample( E, dir, rng );
  }
  else {
    return _rxn_fission_iter->second.scatter_physics->sample( E, dir, rng );
  }
}

//adjust cross sections by factor over a prescribed energy range for testing purposes
//note that values are adjusted in a pointwise manner including both the lower grid point
//containing emin and the upper grid point for the bin containing emax
void xs::neutronData::adjust( const xs::endf_mt_t mt, const double factor, const double emin,
                              const double emax ) {
  auto& rxn     = _rxn.at(mt);
  auto& xs_data = rxn.xs;

  const auto start = emin == 0.0 ? rxn.ierg : util::binary_search( emin, _egrid ) - rxn.ierg;
  const auto end   = emax == std::numeric_limits<double>::max() ? start + rxn.nerg :
    std::min( start + rxn.nerg,  1 + util::binary_search( emax, _egrid ) - rxn.ierg );
  //const auto  start = rxn.ierg;  //start of the portion of energy grid valid for rxn
  //const auto  end   = start + rxn.nerg - 1;  //end of portion of energy grid

  auto& total_rxn = _rxn.at( xs::endf_mt_t::TOTAL );
  auto& total_xs  = total_rxn.xs;

  //loop over reaction compute delta at each energy grid point
  //adjust cross section and preserve total
  for ( auto i = start ; i < end ; ++i ) {
    const auto xs_provided = xs_data[i];
    const auto xs_adjusted = factor * xs_provided;
    const auto xs_change   = xs_adjusted - xs_provided;

    xs_data[i]            = xs_adjusted;
    total_xs[i+rxn.ierg] += xs_change;
  }

  //unresolved resonances also need to be adjusted for elastic, fission, (n,gamma)
  if ( _urr_exist ) {
    _urr->adjust( mt, factor );
  }
}
