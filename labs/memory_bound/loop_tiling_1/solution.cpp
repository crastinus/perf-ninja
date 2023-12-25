#include "solution.hpp"
#include <algorithm>

const int TILE_SIZE = 16;

bool solution(MatrixOfDoubles &in, MatrixOfDoubles &out) {
  int size = in.size();
  for (int i = 0; i < size; i += TILE_SIZE) {
    for (int j = 0; j < size; j += TILE_SIZE) {
      for (int i1 = i; i1 < std::min(size,i + TILE_SIZE); i1++) {
        for (int j1 = j; j1 < std::min(size, j + TILE_SIZE); j1++) {
          out[i1][j1] = in[j1][i1];
        }
      }
    }
  }
  return out[0][size - 1];
}
