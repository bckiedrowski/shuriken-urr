#ifndef _XS_YIELD_HEADER_
#define _XS_YIELD_HEADER_

#include <vector>

#include "ace.hpp"
#include "interpGrid.hpp"

namespace xs {

  //base yield class
  class yieldData {
     public:
       virtual double yield( const double erg ) const = 0;
  };

  class fixedYield : public yieldData {
    private:
      uint32_t _yield;
    public:
      fixedYield( const uint32_t Y ) : _yield(Y) {};

      double yield( const double erg ) const { return _yield; }
  };

  class polynomialYield : public yieldData {
     // yield described by a polynomial in energy
     private:
       std::vector<double> _c;
     public:
       polynomialYield( const aceData& acefile, const uint32_t loc );

       double yield( const double erg ) const;
  };

  class tabularYield : public yieldData {
    // yield described by tabular data
    private:
      interpGrid _yield;
    public:
      tabularYield( const aceData& acefile, const uint32_t loc ) :
        _yield( interpGrid(acefile, loc) ) {};

      double yield( const double erg ) const { return _yield( erg ); };
  };
}

#endif
