#ifndef _RANDOM_HEADER_
#define _RANDOM_HEADER_

#include <cassert>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace util {

  //parent class to hold a true random number generator and a testing class
  class RandomNumberGenerator {
    protected:
      uint64_t hist_id;

      virtual uint64_t skip_ahead(uint64_t* seed, int64_t* nskip) {
        assert(false);
        return 0; };

    public:
      using ptr = std::shared_ptr<RandomNumberGenerator>;

    virtual const std::string getMode() const = 0;
    virtual uint64_t* getSEED0() {
      assert(false);
      return 0;
    };
    virtual uint64_t* getSEED() {
      assert(false);
      return 0;
    };

    virtual double sample() = 0;
    virtual void set(uint64_t nps) = 0;
    virtual uint64_t history() { return hist_id; }

  };

  using RNG = RandomNumberGenerator;

  //standard random number generation (MCNP5 63-bit generator #5) class
  class Rand : public RandomNumberGenerator {
    private:
      int RN_INDEX = 1;
      uint64_t RN_MULT = 3512401965023503517ULL;
      uint64_t RN_ADD = 0ULL;
      int RN_BITS = 63;
      uint64_t RN_STRIDE = 152917ULL;
      uint64_t RN_SEED0 = 1ULL;
      uint64_t RN_MOD = 1ULL << 63;
      uint64_t RN_MASK = (~0ULL) >> 1;
      uint64_t RN_PERIOD = 1ULL << 61;
      double RN_NORM = 1. / static_cast<double>(1ULL << 63);
      uint64_t RN_SEED = 1ULL;

    public:
      uint64_t* getSEED0() override { return &RN_SEED0; };
      uint64_t* getSEED()  override { return &RN_SEED; };
      const std::string getMode() const override { return "mcnpp5 63-bit gen 5"; }

      // call to return a uniform random number between 0 and 1
      double sample() override;

      void set(uint64_t nps) override;
      uint64_t skip_ahead(uint64_t* seed, int64_t* nskip) override;
  };

  //testing random number returns a defined sequence of numbers
  class TestRandomSequence : public RandomNumberGenerator {
    private:
      int testIndex;
      std::vector<double> loopThrough;

    public:
      using ptr = std::shared_ptr<TestRandomSequence>;

      TestRandomSequence(std::vector<double> loopThrough_in)
          : testIndex(0), loopThrough(loopThrough_in){};
      TestRandomSequence() : testIndex(0){};

      const std::string getMode() const override { return "Testing"; }

      //set the vector that the testing mode will return
      void setLoopThrough(std::vector<double> loopThrough_in) {
        loopThrough = loopThrough_in;
        testIndex = 0;
      }

      //set the place in the loopThrough vector
      void setIndex(int index) { testIndex = index; }

      //returns next random number in the sequence
      double sample() override;

      //this function resets the index
      void set(uint64_t nps) override;
  };
} // namespace util

#endif
