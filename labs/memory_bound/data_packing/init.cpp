
#include "solution.h"
#include <random>

// 0 < 100
S create_entry(int first_value, int second_value) {
  S entry;

  entry.i = first_value; // 0..100 8bit
  entry.s = static_cast<short>(second_value); // 0..100 8bit
  entry.l = static_cast<long long>(first_value * second_value); // 0..10000 16bit
  entry.d = static_cast<double>(first_value) / maxRandom; // 0.00..1.00
  entry.b = first_value < second_value; // 1/2  1 bit

  return entry;
}

void init(std::array<S, N> &arr) {
  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution(minRandom, maxRandom - 1);

  for (int i = 0; i < N; i++) {
    int random_int1 = distribution(generator);
    int random_int2 = distribution(generator);

    arr[i] = create_entry(random_int1, random_int2);
  }
}
