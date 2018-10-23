#ifndef BLOOMFILTER_HPP
#define BLOOMFILTER_HPP

typedef uint64_t bf_t;

// l label
// z number of bits
// h number of hash functions
// rng random number generator
bf_t createBloomFilter(int l, int z, int h, std::mt19937 rng)
{
  bf_t out = bf_t(0);
  for(int i=0; i<h; i++) out |= 1<<(rng()%z);
  return out;
}

#endif

