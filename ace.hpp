#ifndef _XS_ACE_HEADER_
#define _XS_ACE_HEADER_

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace xs {

  class aceData {
    private:
      double   _awr;
      double   _temp;

      uint32_t    _zaid;
      std::string _date;
      std::string _info;

      std::array<int,    16> _iz;
      std::array<double, 16> _az;
      std::array<int,    16> _nxs;  // numerical integer data
      std::array<int,    32> _jxs;  // locators in xss vector

      std::vector<double> _xss;     // data
    public:
      aceData( const std::string filename ) { read(filename); }

      // note: use Fortran indexing (starting from 1) as accessor arguments
      //       to be consistent with ace format descriptions
      int nxs( uint32_t k ) const { return _nxs[k-1]; };
      int jxs( uint32_t k ) const { return _jxs[k-1]; };
      double xss( uint32_t k ) const { return _xss[k-1]; };

      double awr()  const { return _awr; }
      double temp() const { return _temp; }

      uint32_t zaid() const { return _zaid; };

      bool read( std::string filename );
      bool write( std::string filename );
  };

}

#endif
