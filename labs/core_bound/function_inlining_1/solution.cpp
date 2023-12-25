
#include "solution.h"
#include <algorithm>
#include <stdlib.h>
#include <iostream>

static int compare(const S& a, const S& b) {

  if (a.key1 < b.key1)
    return -1;

  if (a.key1 > b.key1)
    return 1;

  if (a.key2 < b.key2)
    return -1;

  if (a.key2 > b.key2)
    return 1;

  return 0;
}

void solution(std::array<S, N> &arr) {
  std::sort(arr.begin(), arr.end(), [](const S& a, const S& b) {
      return a.key1 < b.key1 || (a.key1 == b.key1) && (a.key2 < b.key2);
      //return compare(a,b);
  });
  // qsort(arr.data(), arr.size(), sizeof(S), compare);
}
