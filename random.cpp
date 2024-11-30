//63-bit gen #5 from MCNP5 adapted by code published by Forrest Brown in LA-UR-07-7961
#include "random.hpp"

using namespace util;

double TestRandomSequence::sample() {
  double toReturn = loopThrough.at(testIndex);
  testIndex += 1;

  //resets to zero at the end of the loop
  testIndex = testIndex % loopThrough.size();

  return static_cast<double>(toReturn);
}

//initialize in this case just resets the index
void TestRandomSequence::set(uint64_t nps) {
  testIndex = 0;
  hist_id = nps;
}

double Rand::sample() {
  // MCNP random number generator
  RN_SEED = (RN_MULT * RN_SEED) & RN_MASK;
  return static_cast<double>(RN_SEED * RN_NORM);
}

uint64_t Rand::skip_ahead(uint64_t* s, int64_t* n) {
  //  skip ahead n RNs:   RN_SEED*RN_MULT^n mod RN_MOD
  uint64_t seed = *s;
  int64_t nskip = *n;
  while (nskip < 0) {
    nskip += RN_PERIOD; // add period until >0
  }
  nskip = nskip & RN_MASK; // mod RN_MOD
  uint64_t gen = 1, g = RN_MULT, gp, inc = 0, c = RN_ADD, rn;
  // get gen=RN_MULT^n,  in log2(n) ops, not n ops !
  for (; nskip; nskip >>= 1) {
    if (nskip & 1) {
      gen = gen * g & RN_MASK;
      inc = (inc * g + c) & RN_MASK;
    }
    c = g * c + c & RN_MASK;
    g = g * g & RN_MASK;
  }
  rn = (gen * seed + inc) & RN_MASK;

  return static_cast<uint64_t>( rn );
}

void Rand::set(uint64_t nps) {
  // initialize MCNP random number parameters for particle "nps"
  //
  //     * generate a new particle seed from the base seed
  //       & particle index
  //     * set the RN count to zero
  hist_id = nps;
  int64_t nskp = nps * RN_STRIDE;
  RN_SEED = skip_ahead(&RN_SEED0, &nskp);
}
