#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "binarySearch.hpp"
#include "expIntegral.hpp"

#include "urr.hpp"

xs::urrData::urrData( const aceData& acefile ) {
  const auto loc_unr = static_cast<int32_t>( acefile.jxs(23) );

  _nerg           = static_cast<int32_t>(  acefile.xss( loc_unr   ));
  _nprob          = static_cast<int32_t>(  acefile.xss( loc_unr+1 ));
  _interp         = static_cast<interp_t>( acefile.xss( loc_unr+2 ));
  _inelastic_flag = static_cast<int32_t>(  acefile.xss( loc_unr+3 ));
  _other_abs_flag = static_cast<int32_t>(  acefile.xss( loc_unr+4 ));
  _factors_flag   = static_cast<int32_t>(  acefile.xss( loc_unr+5 ));

  assert( _interp == interp_t::LINLIN || _interp == interp_t::LOGLOG );

  _egrid.resize( _nerg );
  _cdf.resize( _nerg );
  _total.resize( _nerg );
  _elastic.resize( _nerg );
  _fission.resize( _nerg );
  _ngamma.resize( _nerg );
  _heating.resize( _nerg );

  _pdf.resize( _nerg );

  const auto loc_egrid  = loc_unr + 6;
  const auto loc_data   = loc_egrid + _nerg;
  const auto size_block = 6 * _nprob;
  for ( auto g = 0 ; g < _nerg ; ++g ) {
    _egrid[g] = acefile.xss( loc_egrid+g );

    //insert zero first element in cdf table to make this work with binary search
    //routine that returns the lower bound.
    //note that ace file only allows discrete values for the probability table
    //if this is changed in the future, this code would need to be rewritten
    _cdf[g].resize( _nprob+1 );
    _cdf[g][0] = 0.0;

    _total[g].resize( _nprob );
    _elastic[g].resize( _nprob );
    _fission[g].resize( _nprob );
    _ngamma[g].resize( _nprob );
    _heating[g].resize( _nprob );
    for ( auto j = 0 ; j < _nprob ; ++j ) {
      const auto k = loc_data + g*size_block + j;
      _cdf[g][j+1]   = acefile.xss( k            ); //note offset on cdf grid
      _total[g][j]   = acefile.xss( k + _nprob   );
      _elastic[g][j] = acefile.xss( k + _nprob*2 );
      _fission[g][j] = acefile.xss( k + _nprob*3 );
      _ngamma[g][j]  = acefile.xss( k + _nprob*4 );
      _heating[g][j] = acefile.xss( k + _nprob*5 );
    }

    //pdf grid (note off by one issue because of cdf search)
    _pdf[g].resize( _nprob );
    for ( auto j = 0 ; j < _nprob ; ++j ) {
      _pdf[g][j] = _cdf[g][j+1] - _cdf[g][j];
    }
  }

  //get information mapping global xs energy onto the unresolved region
  const auto loc_esg = acefile.jxs(1);  //location of energy grid and primary cross sections
  const auto nes     = acefile.nxs(3);  //number of energies on overall energy grid

  //read in overall energy grid and find portion bounding unresolved range
  std::vector< double > overall_egrid( nes );
  for ( auto i = 0 ; i < nes ; ++i ) {
    overall_egrid[i] = acefile.xss( loc_esg + i );
  }
  const auto start_urr = util::binary_search( _egrid[0], overall_egrid );
  const auto end_urr   = util::binary_search( _egrid[_egrid.size()-1], overall_egrid) + 1;
  const auto num_erg_urr = end_urr - start_urr + 1;

  background_xs.egrid.resize( num_erg_urr );
  for ( auto i = 0 ; i < num_erg_urr ; ++i ) {
    background_xs.egrid[i] = overall_egrid[ start_urr + i ];
  }

  //populate background cross sections (only strictly needed for sampling if
  //factors are being used as opposed to values)
  background_xs.elastic.resize( num_erg_urr, 0.0 );
  background_xs.fission.resize( num_erg_urr, 0.0 );
  background_xs.ngamma.resize(  num_erg_urr, 0.0 );
  background_xs.heating.resize( num_erg_urr, 0.0 );

  //elastic and heating values are located in main xs array at head of ace file
  //other reactions are located in other parts of the ace file
  const auto loc_elastic = loc_esg + 3*nes;
  const auto loc_heating = loc_esg + 4*nes;
  for ( auto i = 0 ; i < num_erg_urr ; ++i ) {
    background_xs.elastic[i] = acefile.xss( loc_elastic + start_urr + i - 1 );
    background_xs.heating[i] = acefile.xss( loc_heating + start_urr + i - 1 );
  }

  //search for fission (either MT = 18 or 19, depending on how it is specified)
  //check total fission (MT=18) first
  //note: tabulate_background function handles thresholds and is an issue because
  //      tables may only be given over a specific energy range that may or may not
  //      coincide with all or part of the unresolved range
  auto found_fission = tabulate_background( acefile, start_urr, num_erg_urr,
                                            endf_mt_t::FISSION, background_xs.fission );
  if ( !found_fission ) {
    //if not there, try first-chance fission (MT=19)
    tabulate_background( acefile, start_urr, num_erg_urr,
                         endf_mt_t::Z_F, background_xs.fission );
  }
  //search for (n,gamma) MT = 102
  tabulate_background( acefile, start_urr, num_erg_urr,
                       endf_mt_t::Z_GAMMA, background_xs.ngamma );

  //inelastic cross section
  //to simplify things, this will be precomputed based on MT=51 through 91 regardless of flag
  //exception is a value of -1, which means no inelastic contribution to unresolved range
  //also include MT=5 with inelastic for balance in case (most often this is above range)
  background_xs.inelastic.resize( num_erg_urr, 0.0 );
  if ( _inelastic_flag >= 0 ) {
    //loop until inelastic level not found
    bool more = true;
    for ( auto mt = 51 ; mt < 91 && more ; ++mt ) {
      more = tabulate_background( acefile, start_urr, num_erg_urr,
                                  static_cast<endf_mt_t>(mt), background_xs.inelastic );
    }
    tabulate_background( acefile, start_urr, num_erg_urr,
                         endf_mt_t::Z_N_CONTINUUM, background_xs.inelastic );
  }
  tabulate_background( acefile, start_urr, num_erg_urr,
                       endf_mt_t::Z_ANYTHING, background_xs.inelastic );

  //other absorption reactions (MT = 103 through 117)
  //as with inelastic, if >= 0, precompute by summing over all MT numbers
  background_xs.other_abs.resize( num_erg_urr, 0.0 );
  if ( _other_abs_flag >= 0 ) {
    for ( auto mt = 103 ; mt < 118 ; ++mt ) {
      tabulate_background( acefile, start_urr, num_erg_urr,
                           static_cast<endf_mt_t>(mt), background_xs.other_abs );
    }
  }

}

bool xs::urrData::tabulate_background( const aceData& acefile, const int32_t start_urr, const int32_t num_erg_urr,
                                       const endf_mt_t target_mt, std::vector<double>& xs ) const {
  const auto loc_mtr  = acefile.jxs(3);  //location of MT block
  const auto loc_lsig = acefile.jxs(6);  //location of cross scection locators
  const auto loc_sig  = acefile.jxs(7);  //location of cross section tables
  const auto nmt      = acefile.nxs(4);  //number of reactions excluding elastic

  auto loca = 0;
  for ( auto i = 0 ; i < nmt ; ++i ) {
    const auto mt = static_cast<endf_mt_t>( acefile.xss( loc_mtr + i ) );
    if ( mt == target_mt ) {
      loca = static_cast<int32_t>( acefile.xss( loc_lsig + i ) );
      break;
    }
  }
  //loca must be greater than zero if reaction is found
  if ( loca > 0 ) {
    //start_urr = index at start of energy array corresponding to lower bound of urr region
    //start_rxn = index at start of energy array corresponding to threshold of reaction
    const auto start_rxn = static_cast<int32_t>( acefile.xss(loc_sig + loca - 1) ) - 1;
    const auto num_rxn   = static_cast<int32_t>( acefile.xss(loc_sig + loca ) );

    //offsets provide location for storing and accessing the urr and rxn data respectively
    //note this is complicated by the fact the reaction theshold could be below, within,
    //or above the unresolved range, and also the table could (in principle) stop anywhere.
    //this logic handles all the possible cases
    const auto offset_urr = std::max( 0, start_rxn - start_urr );
    const auto offset_rxn = std::max( 0, start_urr - start_rxn );

    const auto loc_data = loc_sig + loca + 1;
    for ( auto i = 0 ; i + offset_urr < num_erg_urr && i + offset_rxn < num_rxn ; ++i ) {
      xs[ i + offset_urr ] += acefile.xss( loc_data + i + offset_rxn );
    }
  }
  return loca > 0;
}

//compute cross sections given a random number
//note: it is up to the host code to track and store random numbers to preserve correlations
//      between particles and data tables at different temperatures
//note: this should only be called when incident energy is within the unresolved range
xs::urrXS xs::urrData::xs( const double erg, const double random_number ) const {
  assert( random_number >= 0.0 && random_number < 1.0 );
  assert( within_urr_range( erg ) );

  urrXS urr_xs;

  //find index (lower bound) on urr energy grid
  const auto ie = util::binary_search( erg, _egrid );
  urr_xs.energy_index = ie;

  //interpolation parameter between energy grids
  //only LINLIN or LOGLOG supported and checked by assertion in constructor
  const auto rtp = _interp == interp_t::LINLIN ?
    ( erg - _egrid[ie] )/( _egrid[ie+1] - _egrid[ie] ) :
    std::log(erg/_egrid[ie]) / std::log(_egrid[ie+1]/_egrid[ie]);

  urr_xs.interp_factor = rtp;

  //search indices on cdf tables based on provided random number
  const auto ip_bot = util::binary_search( random_number, _cdf[ie]   );
  const auto ip_top = util::binary_search( random_number, _cdf[ie+1] );

  urr_xs.left_xs_index  = ip_bot;
  urr_xs.right_xs_index = ip_top;

  urr_xs.left_elastic_xs  = _elastic[ie][ip_bot];
  urr_xs.left_fission_xs  = _fission[ie][ip_bot];
  urr_xs.left_ngamma_xs   = _ngamma[ie][ip_bot];

  urr_xs.right_elastic_xs  = _elastic[ie+1][ip_top];
  urr_xs.right_fission_xs  = _fission[ie+1][ip_top];
  urr_xs.right_ngamma_xs   = _ngamma[ie+1][ip_top];

  //note: mcnp appears to implicitly use discrete values for the cross sections/factors
  //      on the unresolved probability tables and does not provide a mechanism in the
  //      ace file to change this. perhaps in the future, this could be expanded
  //note: the interpolation over the energy grid (LINLIN or LOGLOG) is specified for
  //      what is provided in the table: cross sections or factors. in the latter case
  //      interpolate factors and then apply background xs after.
  if ( _interp == interp_t::LINLIN ) {
    urr_xs.elastic = _elastic[ie][ip_bot] + rtp*( _elastic[ie+1][ip_top] - _elastic[ie][ip_bot] );
    urr_xs.fission = _fission[ie][ip_bot] + rtp*( _fission[ie+1][ip_top] - _fission[ie][ip_bot] );
    urr_xs.ngamma  = _ngamma[ie][ip_bot]  + rtp*( _ngamma[ie+1][ip_top]  - _ngamma[ie][ip_bot]  );
    urr_xs.heating = _heating[ie][ip_bot] + rtp*( _heating[ie+1][ip_top] - _heating[ie][ip_bot] );
  }
  else {
    urr_xs.elastic = _elastic[ie][ip_bot] > 0.0 ?
      _elastic[ie][ip_bot] * std::pow( _elastic[ie+1][ip_top] / _elastic[ie][ip_bot], rtp ) : 0.0;
    urr_xs.fission = _fission[ie][ip_bot] > 0.0 ?
      _fission[ie][ip_bot] * std::pow( _fission[ie+1][ip_top] / _fission[ie][ip_bot], rtp ) : 0.0;
    urr_xs.ngamma  = _ngamma[ie][ip_bot] > 0.0  ?
      _ngamma[ie][ip_bot]  * std::pow( _ngamma[ie+1][ip_top]  / _ngamma[ie][ip_bot],  rtp ) : 0.0;
    urr_xs.heating = _heating[ie][ip_bot] > 0.0 ?
      _heating[ie][ip_bot] * std::pow( _heating[ie+1][ip_top] / _heating[ie][ip_bot], rtp ) : 0.0;
  }

  //interpolation parameter on portion of global/background xs energy grid in urr range
  //linear interpolation of the background cross sections is assumed
  const auto ib  = util::binary_search( erg, background_xs.egrid );
  const auto rtb = ( erg - background_xs.egrid[ib] )/( background_xs.egrid[ib+1] - background_xs.egrid[ib] );
  if ( _factors_flag > 0 ) {
    //data on table and referenced by xs_ptr are factors and need to multiply background xs
    urr_xs.elastic *= background_xs.elastic[ib] + rtb*( background_xs.elastic[ib+1] - background_xs.elastic[ib] );
    urr_xs.fission *= background_xs.fission[ib] + rtb*( background_xs.fission[ib+1] - background_xs.fission[ib] );
    urr_xs.ngamma  *= background_xs.ngamma[ib]  + rtb*( background_xs.ngamma[ib+1]  - background_xs.ngamma[ib]  );
    urr_xs.heating *= background_xs.heating[ib] + rtb*( background_xs.heating[ib+1] - background_xs.heating[ib] );

    urr_xs.left_elastic_xs *= background_xs.elastic[ib] + rtb*( background_xs.elastic[ib+1] - background_xs.elastic[ib] );
    urr_xs.left_fission_xs *= background_xs.fission[ib] + rtb*( background_xs.fission[ib+1] - background_xs.fission[ib] );
    urr_xs.left_ngamma_xs  *= background_xs.ngamma[ib]  + rtb*( background_xs.ngamma[ib+1]  - background_xs.ngamma[ib]  );

    urr_xs.right_elastic_xs *= background_xs.elastic[ib] + rtb*( background_xs.elastic[ib+1] - background_xs.elastic[ib] );
    urr_xs.right_fission_xs *= background_xs.fission[ib] + rtb*( background_xs.fission[ib+1] - background_xs.fission[ib] );
    urr_xs.right_ngamma_xs  *= background_xs.ngamma[ib]  + rtb*( background_xs.ngamma[ib+1]  - background_xs.ngamma[ib]  );
  }
  //compute total cross section by using sampled values and background data
  const auto inelastic = background_xs.inelastic[ib] + rtb*( background_xs.inelastic[ib+1] - background_xs.inelastic[ib] );
  const auto other_abs = background_xs.other_abs[ib] + rtb*( background_xs.other_abs[ib+1] - background_xs.other_abs[ib] );
  urr_xs.total = urr_xs.elastic + urr_xs.fission + urr_xs.ngamma + inelastic + other_abs;

  urr_xs.other = inelastic + other_abs;

  return urr_xs;
}

std::vector< xs::urrData::union_grid > xs::urrData::unionize() {
  std::vector< union_grid > ugrid;
  ugrid.resize( _nerg - 1 );

  const auto u_nprob = 2*_nprob -1;

  for ( auto g = 0 ; g < _nerg-1 ; ++g ) {
    ugrid[g].i1.resize( u_nprob );
    ugrid[g].i2.resize( u_nprob );
    ugrid[g].pdf.resize( u_nprob );
    ugrid[g].cdf.resize( u_nprob + 1 );

    ugrid[g].cdf[0] = 0.0;

    uint32_t i1 = 0, i2 = 0;
    for ( auto k = 0 ; k < u_nprob ; ++k ) {
      ugrid[g].i1[k] = i1;
      ugrid[g].i2[k] = i2;
      if ( _cdf[g][i1+1] <= _cdf[g+1][i2+1] ) {
        ugrid[g].cdf[k+1] = _cdf[g][i1+1];
        ++i1;
      }
      else {
        ugrid[g].cdf[k+1] = _cdf[g+1][i2+1];
        ++i2;
      }
      ugrid[g].pdf[k] = ugrid[g].cdf[k+1] - ugrid[g].cdf[k];
    }
  }
  return ugrid;
}

void xs::urrData::analytical_benchmark_1( const double N, const double X ) {
  const auto Nx      = N*X;
  const auto ugrid   = unionize();
  const auto u_nprob = 2*_nprob -1;

  std::vector< std::vector<double> > Dleakage_Dsigma;

  Dleakage_Dsigma.resize( _nerg );
  for ( auto g = 0 ; g < _nerg ; ++g ) {
    Dleakage_Dsigma[g].resize( _nprob, 0.0 );
  }

  double leakage = 0.0;
  for ( auto g = 0 ; g < _nerg-1 ; ++g ) {
    const double pg = ( _egrid[g+1] - _egrid[g] )/( upper_bound() - lower_bound() );

    double leak_g = 0.0;
    for ( auto k = 0 ; k < u_nprob ; ++k ) {
      const auto i1 = ugrid[g].i1[k];
      const auto i2 = ugrid[g].i2[k];
      const auto s1 = Nx * ( _elastic[g][i1]   + _fission[g][i1]   + _ngamma[g][i1] );
      const auto s2 = Nx * ( _elastic[g+1][i2] + _fission[g+1][i2] + _ngamma[g+1][i2] );

      const auto g1 = Nx * _ngamma[g][i1];
      const auto g2 = Nx * _ngamma[g+1][i2];

      const auto e1 = Nx * _elastic[g][i1];
      const auto e2 = Nx * _elastic[g+1][i2];

      const auto fleak = ( std::exp(-s1) - std::exp(-s2) )/( s2 - s1 );
      leak_g += fleak * ugrid[g].pdf[k];

      Dleakage_Dsigma[g][i1]     += s1 * ( fleak - std::exp(-s1) )/( s2 - s1 ) * ugrid[g].pdf[k] * pg;
      Dleakage_Dsigma[g+1][i2]   += s2 * ( std::exp(-s2) - fleak )/( s2 - s1 ) * ugrid[g].pdf[k] * pg;
    }
    leakage += leak_g * pg;
  }
  std::cout << "ref. leakage = " << leakage << '\n' << '\n';

  std::cout << "total xs sensitivity matrix = " << '\n';
  double sum_Dleakage_Dsigma = 0.0;
  for ( auto g = 0 ; g < _nerg ; ++g ) {
    for ( auto k = 0 ; k < _nprob ; ++k ) {
      std::cout << std::scientific << std::setw(12) << std::setprecision(4) << Dleakage_Dsigma[g][k] / leakage;
      sum_Dleakage_Dsigma += Dleakage_Dsigma[g][k];
    }
    std::cout << '\n';
  }
  std::cout << std::fixed << std::setprecision(6) << sum_Dleakage_Dsigma / leakage << '\n';

}
void xs::urrData::analytical_benchmark_2( const double N, const double X ) {
  const auto Nx      = N*X;
  const auto ugrid   = unionize();
  const auto u_nprob = 2*_nprob -1;

  std::vector< std::vector<double> > Dcapture_Dgamma;
  std::vector< std::vector<double> > Dcapture_Dscatter;

  Dcapture_Dgamma.resize( _nerg );
  Dcapture_Dscatter.resize( _nerg );
  for ( auto g = 0 ; g < _nerg ; ++g ) {
    Dcapture_Dgamma[g].resize( _nprob, 0.0 );
    Dcapture_Dscatter[g].resize( _nprob, 0.0 );
  }

  double capture = 0.0;
  for ( auto g = 0 ; g < _nerg-1 ; ++g ) {
    const double pg = ( _egrid[g+1] - _egrid[g] )/( upper_bound() - lower_bound() );

    double cap_g  = 0.0;
    for ( auto k = 0 ; k < u_nprob ; ++k ) {
      const auto i1 = ugrid[g].i1[k];
      const auto i2 = ugrid[g].i2[k];
      const auto s1 = Nx * ( _elastic[g][i1]   + _fission[g][i1]   + _ngamma[g][i1] );
      const auto s2 = Nx * ( _elastic[g+1][i2] + _fission[g+1][i2] + _ngamma[g+1][i2] );

      const auto g1 = Nx * _ngamma[g][i1];
      const auto g2 = Nx * _ngamma[g+1][i2];

      const auto e1 = Nx * _elastic[g][i1];
      const auto e2 = Nx * _elastic[g+1][i2];

      const auto hcap = ( g1 - g2 ) * ( s1 - s2 + std::exp(-s1) - std::exp(-s2) )
                      + ( g2*s1 - g1*s2 ) * ( std::log(s1/s2) + util::expint(1,s1) - util::expint(1,s2) );
      const auto fcap = hcap / std::pow( s1 - s2, 2 );
      cap_g += fcap * ugrid[g].pdf[k];

      const auto dhcap1 = s1 - s2 + std::exp(-s1) - std::exp(-s2) + ( g1 - g2 ) * ( 1.0 - std::exp(-s1) )
                        + ( g2 - s2 ) * (  std::log(s1/s2) + util::expint(1,s1) - util::expint(1,s2) )
                        + ( g2*s1 - g1*s2 ) * ( 1.0 - std::exp(-s1) )/s1;
      const auto dhcap2 = s2 - s1 + std::exp(-s2) - std::exp(-s1) + ( g2 - g1 ) * ( 1.0 - std::exp(-s2) )
                        + ( s1 - g1 ) * (  std::log(s1/s2) + util::expint(1,s1) - util::expint(1,s2) )
                        + ( g1*s2 - g2*s1 ) * ( 1.0 - std::exp(-s2) )/s2;

      const auto dhsct1 = ( g1 - g2 ) * ( 1.0 - std::exp(-s1) ) 
                        + g2 * (  std::log(s1/s2) + util::expint(1,s1) - util::expint(1,s2) )
                        + ( g2*s1 - g1*s2 ) * ( 1.0 - std::exp(-s1) )/s1;
      const auto dhsct2 = ( g2 - g1 ) * ( 1.0 - std::exp(-s2) )
                        - g1 * (  std::log(s1/s2) + util::expint(1,s1) - util::expint(1,s2) )
                        + ( g1*s2 - g2*s1 ) * ( 1.0 - std::exp(-s2) )/s2;

      Dcapture_Dgamma[g][i1]     += g1 * ( dhcap1 - 2.0*( s1 - s2 )*fcap ) / std::pow( s1 - s2, 2 ) * ugrid[g].pdf[k] * pg;
      Dcapture_Dgamma[g+1][i2]   += g2 * ( dhcap2 + 2.0*( s1 - s2 )*fcap ) / std::pow( s1 - s2, 2 ) * ugrid[g].pdf[k] * pg;

      Dcapture_Dscatter[g][i1]   += e1 * ( dhsct1 - 2.0*( s1 - s2 )*fcap ) / std::pow( s1 - s2, 2 ) * ugrid[g].pdf[k] * pg;
      Dcapture_Dscatter[g+1][i2] += e2 * ( dhsct2 + 2.0*( s1 - s2 )*fcap ) / std::pow( s1 - s2, 2 ) * ugrid[g].pdf[k] * pg;
    }
    capture += cap_g  * pg;
  }
  std::cout << "ref. capture = " << capture << '\n' << '\n';

  std::cout << "capture xs sensitivity matrix = " << '\n';
  double sum_Dcapture_Dgamma = 0.0;
  for ( auto g = 0 ; g < _nerg ; ++g ) {
    for ( auto k = 0 ; k < _nprob ; ++k ) {
      std::cout << std::scientific << std::setw(12) << std::setprecision(4) << Dcapture_Dgamma[g][k] / capture;
      sum_Dcapture_Dgamma += Dcapture_Dgamma[g][k];
    }
    std::cout << '\n';
  }
  std::cout << std::fixed << std::setprecision(6) << sum_Dcapture_Dgamma / capture << '\n';

  std::cout << "scatter xs sensitivity matrix = " << '\n';
  double sum_Dcapture_Dscatter = 0.0;
  for ( auto g = 0 ; g < _nerg ; ++g ) {
    for ( auto k = 0 ; k < _nprob ; ++k ) {
      std::cout << std::scientific << std::setw(12) << std::setprecision(4) << Dcapture_Dscatter[g][k] / capture;
      sum_Dcapture_Dscatter += Dcapture_Dscatter[g][k];
    }
    std::cout << '\n';
  }
  std::cout << std::fixed << std::setprecision(6) << sum_Dcapture_Dscatter / capture << '\n';
}

void xs::urrData::analytical_benchmark_3( const double N, const double X ) {
  const auto g = 1;
  const auto ugrid   = unionize();
  const auto u_nprob = 2*_nprob -1;

  double leakage = 0.0;

  std::vector< double > Dleakage_Dscatter1( _nprob, 0.0 );
  std::vector< double > Dleakage_Dscatter2( _nprob, 0.0 );

  //std::vector< std::vector<double> > f, df_a, df_b;
  //f.resize( u_nprob ); df_a.resize( u_nprob ); df_b.resize( u_nprob ); 
  for ( auto l = 0 ; l < u_nprob ; ++l ) {
    //f[l].resize( u_nprob, 0.0 ); df_a[l].resize( u_nprob, 0.0 ); df_b[l].resize( u_nprob, 0.0 );

    const auto i1 = ugrid[g].i1[l];
    const auto i2 = ugrid[g].i2[l];
    for ( auto m = 0 ; m < u_nprob ; ++m ) {
      const auto j1 = ugrid[g].i1[m];
      const auto j2 = ugrid[g].i2[m];

      //elastic xs sampled on the upper edge post-collision
      const auto e2f = N * _elastic[g+1][i2];

      //elastic xs sampled on the lower edge postcollision
      const auto e1f = N * _elastic[g][i1];

      //total and elastic cross sections on upper edge sampled pre-collision
      const auto t2i = N * ( _elastic[g+1][j2] + _fission[g+1][j2] + _ngamma[g+1][j2] );
      const auto e2i = N * _elastic[g+1][j2];

      //density times total cross section on lower edge sampled post-collision
      const auto a  = N * ( _elastic[g][i1] + _fission[g][i1] + _ngamma[g][i1] );

      //density times difference in total cross sections on upper edge sampled post-collision minus pre-collision
      const auto b  = N * ( _elastic[g+1][i2] + _fission[g+1][i2] + _ngamma[g+1][i2] )
                    - N * ( _elastic[g+1][j2] + _fission[g+1][j2] + _ngamma[g+1][j2] );

      constexpr auto euler = 0.577215664901533;
      const auto c = ( 1.0 + (b-a)*X ) * std::exp(-(b-a)*X) * util::expint(1,a*X) 
                    + std::log( a/std::abs(b) ) + std::expint(-b*X); 
      const auto h = b == 0 
                   ? -(1.0 - a*X) * std::exp(a*X) * util::expint(1,a*X) - std::log(a*X) - euler  
                   : 1.0 - std::exp(-b*X) + a/(b-a)*c;      
      const auto f = h/(b-a);

      const auto dc_a = X*X*(b-a) * std::exp(-(b-a)*X) * util::expint( 1, a*X )
                      + ( 1.0 - ( 1.0 + (b-a)*X )*std::exp(-b*X) )/a;
      const auto dh_a = b == 0
                      ? -X*( 1 - a * X * std::exp(a*X) * util::expint( 1, a*X ) )
                      : b/std::pow(b-a,2) * c + a/(b-a) * dc_a;
      const auto df_a = h/std::pow( b-a, 2 ) + dh_a/(b-a);


      const auto dc_b = X*X*(a-b) * std::exp(-(b-a)*X) * util::expint( 1, a*X ) - ( 1.0 - std::exp(-b*X) )/b;
      const auto dh_b = b == 0 ? 0.0 : X*std::exp(-b*X) - a*c/std::pow(b-a,2) + a*dc_b/(b-a);
      const auto df_b = dh_b/(b-a) - h/std::pow(b-a,2);

      leakage += ugrid[g].pdf[l] * ugrid[g].pdf[m] * e2i * std::exp(-t2i*X) * f;

      Dleakage_Dscatter1[i1] += ugrid[g].pdf[l] * ugrid[g].pdf[m] * e1f * e2i * std::exp(-t2i*X) * df_a;

      Dleakage_Dscatter2[i2] += ugrid[g].pdf[l] * ugrid[g].pdf[m] * e2f * e2i * std::exp(-t2i*X) * df_b;
      Dleakage_Dscatter2[j2] += ugrid[g].pdf[l] * ugrid[g].pdf[m] * e2i       * std::exp(-t2i*X)
                              * ( ( 1.0 - e2i * X ) * f - e2i * df_b );
    }
  }
  std::cout << " reference once-collided leakage probability: " << leakage << '\n';

  std::cout << " 54 keV sensitivity coefficients: " << '\n';
  double sum_Dleakage_Dscatter1 = 0.0;
  for ( auto k = 0 ; k < _nprob ; ++k ) {
    std::cout << std::scientific << std::setw(12) << std::setprecision(4) << Dleakage_Dscatter1[k] / leakage;
    sum_Dleakage_Dscatter1 += Dleakage_Dscatter1[k];
  }
  std::cout << "\n sum = " << std::fixed << std::setprecision(6) << sum_Dleakage_Dscatter1 / leakage << '\n';
  std::cout << " 59 keV sensitivity coefficients: " << '\n';
  double sum_Dleakage_Dscatter2 = 0.0;
  for ( auto k = 0 ; k < _nprob ; ++k ) {
    std::cout << std::scientific << std::setw(12) << std::setprecision(4) << Dleakage_Dscatter2[k] / leakage;
    sum_Dleakage_Dscatter2 += Dleakage_Dscatter2[k];
  }
  std::cout << "\n sum = " << std::fixed << std::setprecision(6) << sum_Dleakage_Dscatter2 / leakage << '\n';

}

//note this is only really useful if values given as cross sections and not factors
std::pair<double,double> xs::urrData::table_value( const endf_mt_t rxn, const uint32_t energy_index,
                                                   const uint32_t left_index, const uint32_t right_index ) const {

  const auto ie = energy_index;
  const auto il = left_index;
  const auto ir = right_index;

  switch( rxn ) {
  case endf_mt_t::TOTAL:
    return std::make_pair( _elastic[ie  ][il] + _fission[ie  ][il] + _ngamma[ie  ][il],
                           _elastic[ie+1][ir] + _fission[ie+1][ir] + _ngamma[ie+1][ir]  );
  case endf_mt_t::ELASTIC:
    return std::make_pair( _elastic[ie][il], _elastic[ie+1][ir] );
  case endf_mt_t::FISSION:
    return std::make_pair( _fission[ie][il], _fission[ie+1][ir] );
  case endf_mt_t::Z_GAMMA:
    return std::make_pair( _ngamma[ie][il],  _ngamma[ie+1][ir] );
  default:
    return std::make_pair( 0.0, 0.0 );
  }
}

void xs::urrData::adjust( const xs::endf_mt_t mt, const double factor ) {
  //adjust cross sections if table values are cross sections and not as factors
  if ( _factors_flag == 0 ) {
    switch( mt ) {
    case xs::endf_mt_t::ELASTIC:
      for ( auto g = 0 ; g < _nerg ; ++g ) {
        for ( auto k = 0 ; k < _nprob ; ++k ) {
          _elastic[g][k] *= factor;
        }
      }
      break;
    case xs::endf_mt_t::FISSION:
      for ( auto g = 0 ; g < _nerg ; ++g ) {
        for ( auto k = 0 ; k < _nprob ; ++k ) {
          _fission[g][k] *= factor;
        }
      }
      break;
    case xs::endf_mt_t::Z_GAMMA:
      for ( auto g = 0 ; g < _nerg ; ++g ) {
        for ( auto k = 0 ; k < _nprob ; ++k ) {
          _ngamma[g][k] *= factor;
        }
      }
      break;
    default:
      //do nothing
      break;
    }
  }
}
