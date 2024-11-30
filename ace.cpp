#include <fstream>

#include "constants.hpp"

#include "ace.hpp"

bool xs::aceData::read( std::string filename ) {

  std::ifstream ace_file;
  ace_file.open( filename );

  std::string line;

  //skip first three lines (assume 2.0.0 format)
  std::getline( ace_file, line );
  std::getline( ace_file, line );
  std::getline( ace_file, line );

  //read and extract header information
  std::getline( ace_file, line );

  _zaid = std::stoi( line.substr(  0, 10 ) );
  _awr  = std::stod( line.substr( 10, 12 ) );
  _temp = std::stod( line.substr( 22, 12 ) ) / constants::k_boltzmann;  //convert MeV to K
  _date = line.substr( 34, 11 );

  //read and extract information line
  std::getline( ace_file, line );
  _info = line.substr( 0, 80 );

  //read next four lines containing pairs of iz, az data
  { auto k = 0;
  for ( auto i = 0; i < 4 ; ++i ) {
    std::getline( ace_file, line );
    auto l = 0;
    for ( auto j = 0 ; j < 4 ; ++j ) {
      _iz[k] = std::stoi( line.substr( l,   7  ) );
      _az[k] = std::stod( line.substr( l+7, 11 ) );
      l += 18;
      ++k;
    }
  }}

  //read next two lines for nxs array
  { auto k = 0;
  for ( auto i = 0; i < 2 ; ++i ) {
    std::getline( ace_file, line );
    auto l = 0;
    for ( auto j = 0 ; j < 8 ; ++j ) {
      _nxs[k] = std::stoi( line.substr( l, 9 ) );
      l += 9;
      ++k;
    }
  }}

  //read next four lines for jxs array
  { auto k = 0;
  for ( auto i = 0; i < 4 ; ++i ) {
    std::getline( ace_file, line );
    auto l = 0;
    for ( auto j = 0 ; j < 8 ; ++j ) {
      _jxs[k] = std::stoi( line.substr( l, 9 ) );
      l += 9;
      ++k;
    }
  }}

  //read in xss array
  {
  const auto nxss = _nxs[0]; //size = _nxs[0]
  _xss.resize( nxss );

  auto k = 0;
  const auto nlines = 1 + nxss / 4;
  for ( auto i = 0; i < nlines ; ++i ) {
    std::getline( ace_file, line );
    auto l = 0;
    // compute number of entries on the line (normally 4 except on last line)
    const auto nentries = std::min( 4, nxss-i*4 );
    for ( auto j = 0 ; j < nentries ; ++j ) {
      _xss[k] = std::stod( line.substr( l, 20 ) );
      l += 20;
      ++k;
    }
  }}

  ace_file.close();
  return true;  //success
}
