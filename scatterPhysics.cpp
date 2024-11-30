#include <cmath>
#include <iostream>
#include <utility>

#include "binarySearch.hpp"
#include "rotate.hpp"

#include "scatterPhysics.hpp"

//-----------------------------------------------------------------------------
xs::scatterPhysics::energyGridInterp::energyGridInterp( const aceData& acefile, const uint32_t loc ) {
  const auto num_interp = static_cast<int32_t>( acefile.xss(loc) );
  const auto num_ierg   = static_cast<int32_t>( acefile.xss(loc + 2*num_interp + 1) );

  //set up coarse interpolation grid
  if ( num_interp == 0 ) {
    //special case of a single bin with LINLIN
    interp_bounds.resize(1);
    interp_bounds[0] = num_ierg-1;

    interp_types.resize(1);
    interp_types[0] = interp_t::LINLIN;
  }
  else {
    interp_bounds.resize(num_interp);
    interp_types.resize(num_interp);
    for ( auto i = 0 ; i < num_interp ; ++i ) {
      interp_bounds[i] = static_cast<int32_t>( acefile.xss(loc + i + 1) - 1 );
      interp_types[i]  = static_cast<interp_t>( acefile.xss(loc + num_interp + i + 1) );
      assert( interp_types[i] == interp_t::HISTOGRAM || interp_types[i] == interp_t::LINLIN );
    }
  }

  //set up fine incident energy grid
  const auto loc_egrid = loc + 2*num_interp + 2;
  egrid.resize( num_ierg );
  for ( auto i = 0 ; i < num_ierg ; ++i ) {
    egrid[i] = acefile.xss( loc_egrid + i );
  }
}

std::pair< uint32_t, double > xs::scatterPhysics::energyGridInterp::search( const double erg ) const {
  //check out of bounds, return value endpoint
  if ( erg <  egrid[0] ) { return std::make_pair( 0, 0.0 ); }
  if ( erg >= egrid[egrid.size()-1] ) { return std::make_pair( egrid.size()-2, 1.0 ); }

  //binary search on fine grid of values
  const auto ibot = util::binary_search( erg, egrid );
  const auto itop = ibot+1;

  //grab interpolation scheme from coarse grid with linear search since there
  //are usually very few (usually one) grid point
  uint32_t k = 0;
  for ( k = 0 ; k < interp_bounds.size() ; ++k ) {
    if ( itop <= interp_bounds[k] ) break;
  }
  const auto interp = interp_types[k];
  const auto r = interp == interp_t::LINLIN ? (erg - egrid[ibot])/(egrid[itop] - egrid[ibot]) : 0.0;
  //return std::tie( ibot, r );
  return std::pair< uint32_t, double >( ibot, r );
}

//-----------------------------------------------------------------------------
std::pair<double,double> xs::scatterPhysics::aceLaw::com_to_lab(
  const double erg_in, const double erg_cm, const double mu_cm,  const double awr ) const {

  const auto erg_lab = erg_cm + ( erg_in + 2.0*mu_cm*(awr+1.0)*std::sqrt( erg_in*erg_cm ) )/std::pow(awr+1.0,2);
  const auto mu_lab  = mu_cm * std::sqrt( erg_cm / erg_lab ) + std::sqrt( erg_in / erg_lab )/(awr+1.0);

  //return std::tie( erg_lab, mu_lab );
  return std::pair<double,double>( erg_lab, mu_lab );
}

//-----------------------------------------------------------------------------
xs::scatterPhysics::freeGasLaw::freeGasLaw(
  const aceData& acefile, std::shared_ptr<angularDistributionTable> ang ) {

  ang_dist = ang;
  A    = acefile.awr();
  temp = acefile.temp();
}

std::pair<double,util::point> xs::scatterPhysics::freeGasLaw::sample(
  const double erg, const util::point& dir, util::RNG::ptr rng ) {

  if ( A > 1.0 && erg > 4000.0 * constants::k_boltzmann * temp ) {
    //stationary target physics
    const double mu_cm = ang_dist->sample( erg, rng );

    const double alpha = std::pow( (A-1)/(A+1), 2 );
    const double E  = 0.5*erg*( ( 1 - alpha )*mu_cm + 1 + alpha );
    const double mu = ( 1 + A*mu_cm )/std::sqrt( 1 + 2*A*mu_cm + A*A );
    const double azi = constants::two_pi*rng->sample();

    return std::make_pair( E, util::rotate( dir, mu, azi ) );
  }
  else {
    //thermal motion
    //TODO: no doppler broadening resonance correction included
    const auto speed_n_lab = std::sqrt( erg );     //proportional to particle speed (factor of 1/2 divides away)
    const auto vel_n_lab   = dir * speed_n_lab;

    //sample velocity of the target nucleus
    const auto beta_v_n = sqrt( erg*A/(constants::k_boltzmann * temp) );
    const auto p_pdf1   = 2.0/(2.0 + constants::sqrt_pi*beta_v_n);

    double r1, p_accept;
    double mu_target, beta_v_target_sq;
    do {
      r1 = rng->sample();
      const auto r2 = rng->sample();
      if ( r2 < p_pdf1 ) {
        const auto r3 = rng->sample();
        const auto r4 = rng->sample();
        beta_v_target_sq = -std::log(r3*r4);
      }
      else {
        const auto r3 = rng->sample();
        const auto r4 = rng->sample();
        const auto r5 = rng->sample();
        const auto cs = std::cos(constants::half_pi*r5);
        beta_v_target_sq = -log(r3) - log(r4)*cs*cs;
      }
      const auto beta_v_target = std::sqrt( beta_v_target_sq );

      mu_target = 2.0*rng->sample() - 1.0;
      p_accept = std::sqrt( beta_v_n*beta_v_n + beta_v_target_sq - 2.0*beta_v_n*beta_v_target*mu_target ) /
                            ( beta_v_n + beta_v_target );
    } while ( r1 >= p_accept );
    const auto speed_target = std::sqrt( beta_v_target_sq * constants::k_boltzmann * temp / A  );

    const auto azi_target = constants::two_pi * rng->sample();
    const auto vel_target = util::rotate( dir, mu_target, azi_target ) * speed_target;

    const auto vel_com = ( vel_n_lab + vel_target * A )/(1.0 + A);
    const auto vel_n_com = vel_n_lab - vel_com;

    const auto speed_n_com = std::sqrt( util::dot_product( vel_n_com, vel_n_com ) );

    //get center-of-mass scattering cosine and azimuthal angle
    const auto mu_com = ang_dist->sample( erg, rng );
    const auto azi    = constants::two_pi * rng->sample();

    //determine new velocity in center of mass frame and convert back to lab
    const auto dir_n_com = vel_n_com / speed_n_com;
    const auto vel_n_com_out = util::rotate( dir_n_com, mu_com, azi ) * speed_n_com;
    const auto vel_n_lab_out = vel_n_com_out + vel_com;

    //get outgoing lab frame energy and direction vector
    const auto E_lab   = util::dot_product( vel_n_lab_out, vel_n_lab_out );
    const auto dir_lab = vel_n_lab_out / std::sqrt(E_lab);

    return std::make_pair( E_lab, dir_lab );
  }
}

//-----------------------------------------------------------------------------
xs::scatterPhysics::aceLaw3::aceLaw3(
  const aceData& acefile, const uint32_t loc, ref_frame_t frame, std::shared_ptr<angularDistributionTable> ang ) {
  //requires COM frame for discrete level inelastic scattering
  assert( frame == ref_frame_t::COM_FRAME );

  ang_dist = ang;
  A  = acefile.awr();
  t1 = acefile.xss( loc );
  t2 = acefile.xss( loc+1 );
}

std::pair<double,util::point> xs::scatterPhysics::aceLaw3::sample(
  const double erg, const util::point& dir, util::RNG::ptr rng ) {
  const auto mu_cm = ang_dist->sample( erg, rng );
  const auto E_cm  = t2*( erg - t1 );

  const auto [ E, mu ] = com_to_lab( erg, E_cm, mu_cm, A );
  const auto azi = constants::two_pi*rng->sample();
  return std::make_pair( E, util::rotate( dir, mu, azi ) );
}

//-----------------------------------------------------------------------------
xs::scatterPhysics::aceLawTabular::aceLawTabular(
  ace_tabular_t type, const aceData& acefile, const uint32_t loc_dlw, const uint32_t loc_ldat, ref_frame_t frame,
  std::shared_ptr<angularDistributionTable> ang ) :
  law_type(type), inc_egrid( energyGridInterp( acefile, loc_ldat ) ) {

  ref_frame = frame;
  A = acefile.awr();

  const auto nr    = static_cast<uint32_t>( acefile.xss(loc_ldat) );
  const auto nierg = static_cast<uint32_t>( acefile.xss(loc_ldat + 2*nr + 1 ) );

  //setup locators (relative to start of dlw or dlwp block in ace file)
  locators.resize( nierg );
  const auto loc_locators = loc_ldat + 2*nr + nierg + 2;
  for ( auto i = 0 ; i < nierg ; ++i ) {
    locators[i] = loc_dlw + static_cast<uint32_t>( acefile.xss( loc_locators+i ) ) - 1;
  }

  //find number of discrete lines (must be identical for every table)
  //encoded within an integer in the first element of the block (see pg. F-22 of mcnp5 manual, vol III)
  // iint_prime = ( ndiscrete * 10 ) + interpolation flag
  auto iint_prime = static_cast<uint32_t>( acefile.xss( locators[0] ) );
  ndiscrete = iint_prime / 10;

  //construct outgoing energy data
  edata.resize( nierg );
  for ( auto i = 0 ; i < nierg ; ++i ) {
    //interpolation scheme
    iint_prime = static_cast<uint32_t>( acefile.xss( locators[i] ) );
    edata[i].interp = iint_prime % 10 == 1 ? interp_t::HISTOGRAM : interp_t::LINLIN;

    const auto num_points = static_cast<uint32_t>( acefile.xss( locators[i] + 1 ) );
    edata[i].eout.resize( num_points );
    edata[i].pdf.resize(  num_points );
    edata[i].cdf.resize(  num_points );
    for ( auto j = 0 ; j < num_points ; ++j ) {
      edata[i].eout[j] = acefile.xss( locators[i] + 2 + j );
      edata[i].pdf[j]  = acefile.xss( locators[i] + 2 +   num_points + j );
      edata[i].cdf[j]  = acefile.xss( locators[i] + 2 + 2*num_points + j );
    }
  }

  //assign angular information
  switch( law_type ) {
  case ace_tabular_t::LAW4:
    ang_dist  = ang;
    break;
  case ace_tabular_t::LAW44:
    for ( auto i = 0 ; i < nierg ; ++i ) {
      const auto num_points = static_cast<uint32_t>( acefile.xss( locators[i] + 1 ) );
      edata[i].R.resize( num_points );
      edata[i].A.resize( num_points );
      for ( auto j = 0 ; j < num_points ; ++j ) {
        edata[i].R[j] = acefile.xss( locators[i] + 2 + 3*num_points + j );
        edata[i].A[j] = acefile.xss( locators[i] + 2 + 4*num_points + j );
      }
    }
    break;
  case ace_tabular_t::LAW61:
    for ( auto i = 0 ; i < nierg ; ++i ) {
      const auto num_points = static_cast<uint32_t>( acefile.xss( locators[i] + 1 ) );
      edata[i].ang_dist.resize( num_points );
      for ( auto j = 0 ; j < num_points ; ++j ) {
        const auto lc = static_cast<int32_t>( acefile.xss( locators[i] + 2 + 3*num_points + j ) );
        assert( lc >= 0 );
        if ( lc == 0 ) {
          edata[i].ang_dist[j] = std::make_shared<isotropicAngDist>();
        }
        else {
          edata[i].ang_dist[j] = std::make_shared<tabularAngDist>( acefile, loc_dlw + std::abs(lc) - 1 );
        }
      }
    }
    break;
  default:
    //unknown law
    assert(false);
  }
}

std::pair<double,util::point> xs::scatterPhysics::aceLawTabular::sample(
  const double erg, const util::point& dir, util::RNG::ptr rng ) {
  double mu, E;

  //search incident energy grid to find lower index for incident energy table
  //and interpolation fraction with respect to the next table
  const auto [ ie, rtp ] = inc_egrid.search( erg );
  auto le = ie;   //index for incident energy table actually sampled
  auto k  = 0;    //index for outgoing energy grid (lower if lin-lin interpolation)

  //sample whether the emission is in the discrete or continuous portion
  //by checking if random number below the continuous part of the table
  const auto r1 = rng->sample();
  bool within_discrete = false;

  if ( ndiscrete > 0 ) {
    //TODO: could be off by 1 here, inspect a table to verify and remove assertion after fixed
    assert( false );
    const auto c = edata[le].cdf[ndiscrete] + rtp*( edata[le+1].cdf[ndiscrete] - edata[le].cdf[ndiscrete] );
    if ( r1 < c ) {
      //within discrete part of the table, do linear search to find discrete line
      //cannot effectively use binary search since values depend on interpolation fraction
      within_discrete = true;
      while ( r1 < edata[le].cdf[k] + rtp*( edata[le+1].cdf[k] - edata[ie].cdf[k] ) ) { ++k; }
      E = edata[le].eout[k] + rtp*( edata[le+1].eout[k] - edata[le].eout[k] );
    }
  }

  if ( ! within_discrete ) {
    //in continuous part of the table, sample lower or upper table based on interpolation fraction
    //and search that table (increment le from ie to ie + 1 if upper table)
    const auto r2 = rng->sample();
    if ( r2 < rtp ) le++;

    //perform binary search on cdf
    //TODO: check start = ndiscrete is not off by one by looking at example of table with lines
    k = util::binary_search( r1, edata[le].cdf, ndiscrete );
    if ( edata[le].interp == interp_t::HISTOGRAM ) {
      E = edata[le].eout[k] + ( r1 - edata[le].cdf[k] )/edata[le].pdf[k];
    }
    else {
      //LINLIN interpolation
      const auto f = ( edata[le].pdf[k+1] - edata[le].pdf[k] )/( edata[le].eout[k+1] - edata[le].eout[k] );
      if ( std::fabs(f) > 1.0e-9 ) {
        const auto p = edata[le].pdf[k];
        const auto c = edata[le].cdf[k];
        E = edata[le].eout[k] + ( std::sqrt( p*p + 2*f*(r1 - c) ) - p )/f;
      }
      else {
        //pdf is flat (or effectively so) in this range and the formula above is invalid
        //and a histogram treatment is correct or a good approximation
        E = edata[le].eout[k] + ( r1 - edata[le].cdf[k] )/edata[le].pdf[k];
      }
    }
  }

  //sample angular distribution based on law type and what happened above
  if ( law_type == ace_tabular_t::LAW4 ) {
    mu = ang_dist->sample( erg, rng );
  }
  else if ( law_type == ace_tabular_t::LAW44 ) {
    //interpolate get precompound and pdf slopes
    double kalbach_R, kalbach_A;
    if ( within_discrete ) {
      kalbach_R = edata[ie].R[k] + rtp*( edata[ie+1].R[k] - edata[ie].R[k] );
      kalbach_A = edata[ie].A[k] + rtp*( edata[ie+1].A[k] - edata[ie].A[k] );
    }
    else {
      //note: mcnp5 manual is a bit ambiguous (typo?) here about which grid is being interpolated over
      //      it seems most logical to be consistent with what is done above for the continuum
      //      and use the grid actually sampled, not the lower bound for histogrm
      const auto f = edata[le].interp == interp_t::HISTOGRAM
                   ? 0.0 : ( E - edata[le].eout[k] )/( edata[le].eout[k+1] - edata[le].eout[k] );
      kalbach_R = edata[le].R[k] + f*( edata[le].R[k+1] - edata[le].R[k] );
      kalbach_A = edata[le].A[k] + f*( edata[le].A[k+1] - edata[le].A[k] );
    }
    const auto r3 = rng->sample();
    const auto r4 = rng->sample();
    if ( r3 < kalbach_R ) {
      mu = std::log( r4*std::exp(kalbach_A) + (1-r4)*std::exp(-kalbach_A))/kalbach_A;
    }
    else {
      const auto t = ( 2*r4 - 1 )*std::sinh(kalbach_A);
      mu = std::log( t + std::sqrt( t*t + 1 ) )/kalbach_A;
    }
  }
  else if ( law_type == ace_tabular_t::LAW61 ) {
    //note: mcnp5 manual ambiguous how to treat discrete case
    //      decide to take angular distribution on closest incident energy for kth line
    //      this is consistent with the treatment prescribed for lin-lin interpolation
    //      on the outgoing energy grid of the continuum case
    if ( within_discrete ) {
      const auto iu = rtp < 0.5 ? ie : ie+1;
      mu = edata[iu].ang_dist[k]->sample( rng );
    }
    else {
      if ( edata[le].interp == interp_t::HISTOGRAM ) {
        mu = edata[le].ang_dist[k]->sample( rng );
      }
      else {
        //take closest according to cdf for lin-lin
        const auto f  = ( r1 - edata[le].cdf[k] )/( edata[le].cdf[k+1] - edata[le].cdf[k] );
        const auto iu = f < 0.5 ? k : k+1;
        mu = edata[le].ang_dist[iu]->sample( rng );
      }
    }
  }
  else {
    //error (should not hit this)
    assert( false );
  }

  //perform unit-base interpolation on incident energy grid to preserve thresholds
  //if sampled in continuous part of table
  if ( ! within_discrete ) {
    const auto nbot = ndiscrete;

    //outgoing energy arrays of different sizes
    const auto ntop_lower   = edata[ie  ].eout.size() - 1;
    const auto ntop_upper   = edata[ie+1].eout.size() - 1;
    const auto ntop_sampled = edata[le  ].eout.size() - 1;

    const auto Emin = edata[ie].eout[nbot]       + rtp*( edata[ie+1].eout[nbot]       - edata[ie].eout[nbot]       );
    const auto Emax = edata[ie].eout[ntop_lower] + rtp*( edata[ie+1].eout[ntop_upper] - edata[ie].eout[ntop_lower] );
    const auto f = ( E - edata[le].eout[nbot] )/( edata[le].eout[ntop_sampled] - edata[le].eout[nbot] );
    E = Emin + f*( Emax - Emin );
  }

  //transform to lab frame if needed
  const auto azi = constants::two_pi*rng->sample();
  if ( ref_frame == ref_frame_t::COM_FRAME ) {
    const auto [ E_lab, mu_lab ] = com_to_lab( erg, E, mu, A );
    return std::make_pair( E_lab, util::rotate( dir, mu_lab, azi ) );
  }
  return std::make_pair( E, util::rotate( dir, mu, azi ) );
}

//-----------------------------------------------------------------------------
xs::scatterPhysics::aceLaw7::aceLaw7(
  const aceData& acefile, const uint32_t loc, ref_frame_t frame, std::shared_ptr<angularDistributionTable> ang ) :
  interp( interpGrid( acefile, loc ) )  {

  ang_dist  = ang,
  ref_frame = frame;
  A = acefile.awr();

  const auto nr = static_cast<int32_t>( acefile.xss(loc) );
  const auto ne = static_cast<int32_t>( acefile.xss(loc + 1 + 2*nr) );
  restrict_energy = acefile.xss( loc + 2 + 2*nr + 2*ne );
}

std::pair<double,util::point> xs::scatterPhysics::aceLaw7::sample(
  const double erg, const util::point& dir, util::RNG::ptr rng ) {

  double mu = ang_dist->sample( erg, rng );

  double E;
  do {
    double r1, r2, s;
    do {
     r1 = rng->sample();
     r2 = rng->sample();
     s  = r1*r1 + r2*r2;
    } while ( s > 1.0 );
    const auto r3 = rng->sample();
    const auto r4 = rng->sample();
    const auto T  = interp( erg );

    E = -T * ( r1*r1/s*std::log(r3) + std::log(r4) );
  } while ( E > erg - restrict_energy );

  const auto azi = constants::two_pi*rng->sample();
  if ( ref_frame == ref_frame_t::COM_FRAME ) {
    const auto [ E_lab, mu_lab ] = com_to_lab( erg, E, mu, A );
    std::make_pair( E_lab, util::rotate( dir, mu_lab, azi ) );
  }
  return std::make_pair( E, util::rotate( dir, mu, azi ) );
}

//-----------------------------------------------------------------------------
xs::scatterPhysics::aceLaw9::aceLaw9(
  const aceData& acefile, const uint32_t loc, ref_frame_t frame, std::shared_ptr<angularDistributionTable> ang ) :
  interp( interpGrid( acefile, loc ) )  {

  ang_dist  = ang,
  ref_frame = frame;
  A = acefile.awr();

  const auto nr = static_cast<int32_t>( acefile.xss(loc) );
  const auto ne = static_cast<int32_t>( acefile.xss(loc + 1 + 2*nr) );
  restrict_energy = acefile.xss( loc + 2 + 2*nr + 2*ne );
}

std::pair<double,util::point> xs::scatterPhysics::aceLaw9::sample(
  const double erg, const util::point& dir, util::RNG::ptr rng ) {

  double mu = ang_dist->sample( erg, rng );

  double E;
  do {
    const auto r1 = rng->sample();
    const auto r2 = rng->sample();
    const auto T  = interp( erg );

    E = -T * std::log( r1 * r2 );
  } while ( E > erg - restrict_energy );

  const auto azi = constants::two_pi*rng->sample();
  if ( ref_frame == ref_frame_t::COM_FRAME ) {
    const auto [ E_lab, mu_lab ] = com_to_lab( erg, E, mu, A );
    std::make_pair( E_lab, util::rotate( dir, mu_lab, azi ) );
  }
  return std::make_pair( E, util::rotate( dir, mu, azi ) );
}

//-----------------------------------------------------------------------------
xs::scatterPhysics::scatterPhysics( const aceData& acefile, const uint32_t mt_offset ) {
  const auto loc_tyr    = acefile.jxs(5);  //location of yield block (encodes ref. frame)
  const auto loc_land   = acefile.jxs(8);  //location of angular distribution locators
  const auto loc_and    = acefile.jxs(9);  //location of angular distribution block

  _ang_dist = std::make_shared< angularDistributionTable >( acefile, mt_offset );

  //mt_offset = 0 corresponds to elastic, which has a special scattering law
  //angular distribution must be provided in center of mass frame
  //outgoing energy determined from analytical free-gas law
  if ( mt_offset == 0 ) {
    _elastic   = true;
    _ref_frame = ref_frame_t::COM_FRAME;

    lawBlock law_block;

    law_block.id   = 0;
    law_block.data = std::make_shared<freeGasLaw>( acefile, _ang_dist );
    _laws.push_back( std::move(law_block) );
  }
  else {
    //inelastic energy distribution
    //TODO: jxs locations depend on neutron, photon, or delayed emission so cannot use the
    //next two lines of code for that (consider adding a context flag to specify?)
    const auto loc_ldlw = acefile.jxs(10);
    const auto loc_dlw  = acefile.jxs(11);

    //grab reference frame (encoded as sign of reaction yield, see pg. F-16 of mcnp5 manual vol. III)
    const auto ty = static_cast<int32_t>( acefile.xss( loc_tyr + mt_offset - 1 ) );
    assert( ty != 0 );
    _ref_frame = ty < 0 ? ref_frame_t::COM_FRAME : ref_frame_t::LAB_FRAME;
//    std::cout << static_cast<int32_t>( _ref_frame ) << " ";

    create_laws( acefile, loc_ldlw, loc_dlw, mt_offset );
  }
}

std::pair<double,util::point> xs::scatterPhysics::sample(
  const double erg, const util::point& dir, util::RNG::ptr rng ) const {
  //sample the scattering law (if there are multiple) based on prescribed probability
  auto ilaw = 0;
  if ( _laws.size() > 1 ) {
    auto r = rng->sample();
    for ( ; ilaw < _laws.size()-1 ; ++ilaw ) {
      //reduce random number by probability and take once falls below zero
      //note loop index is such that the final law is taken otherwise (avoids roundoff issues)
      r -= _laws[ilaw].prob->operator()( erg );
      if ( r <= 0.0 ) break;
    }
  }
  return _laws[ilaw].data->sample( erg, dir, rng );
}

void xs::scatterPhysics::create_laws( const aceData& acefile, const uint32_t loc_ldlw,
                                      const uint32_t loc_dlw, const uint32_t mt_offset ) {

    const auto locc = static_cast<uint32_t>( acefile.xss( loc_ldlw + mt_offset - 1 ) );
    uint32_t next_law_offset = locc;
    do {
      lawBlock law_block;

      const auto loc_header = loc_dlw + next_law_offset - 1;

      next_law_offset = static_cast<uint32_t>( acefile.xss( loc_header ) );
      if ( next_law_offset > 0 ) {
        //next_law = 0 denotes that this is the last available law for this reaction
        //if this is not the case, then need a probability grid to use that law
        //if the last law is reached in sampling, then this law will be automatically used
        //and the probability table is unnecessary
        law_block.prob = std::make_shared<interpGrid>( acefile, loc_header + 3 );
      }

      law_block.id   = static_cast<uint32_t>( acefile.xss( loc_header + 1 ) );
      const uint32_t loc_data = static_cast<uint32_t>( acefile.xss( loc_header + 2 ) ) + loc_dlw - 1;

      //std::cout << law_block.id << " ";

      switch( law_block.id ) {
      case 3:
        law_block.data = std::make_shared<aceLaw3>( acefile, loc_data, _ref_frame, _ang_dist );
        break;
      case 4:
        law_block.data = std::make_shared<aceLawTabular>(
          ace_tabular_t::LAW4, acefile, loc_dlw, loc_data, _ref_frame, _ang_dist );
        break;
      case 7:
        law_block.data = std::make_shared<aceLaw7>( acefile, loc_data, _ref_frame, _ang_dist );
        break;
      case 9:
        law_block.data = std::make_shared<aceLaw9>( acefile, loc_data, _ref_frame, _ang_dist );
        break;
      case 44:
        law_block.data = std::make_shared<aceLawTabular>(
          ace_tabular_t::LAW44, acefile, loc_dlw, loc_data, _ref_frame );
        break;
      case 61:
        law_block.data = std::make_shared<aceLawTabular>(
          ace_tabular_t::LAW61, acefile, loc_dlw, loc_data, _ref_frame );
        break;
      default:
        assert( false );
      }
      //note: typically the number of laws is very small and the size cannot quickly be
      //      computed ahead of time, so push_back is correct in this instance
      _laws.push_back( std::move(law_block) );
    } while ( next_law_offset > 0 );
}
