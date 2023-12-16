
#include "solution.h"
#include <memory>
#include <string_view>

// Make zero matrix
void zero(Matrix &result) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      result[i][j] = 0;
    }
  }
}

// Make identity matrix
void identity(Matrix &result) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      result[i][j] = 0;
    }
    result[i][i] = 1;
  }
}

void tiling_multiply(Matrix &result, const Matrix &a, const Matrix &b) ;
// Multiply two square matrices
void multiply(Matrix &result, const Matrix &a, const Matrix &b) {
  zero(result);

  //tiling_multiply(result, a, b);

   for (int i = 0; i < N; i++) {
     for (int k = 0; k < N; k++) {
       for (int j = 0; j < N; j++) {
         result[i][j] += a[i][k] * b[k][j];
       }
     }
   }
}

// void tiling_multiply(Matrix &result, const Matrix &a, const Matrix &b) {
//   const int TILE_SIZE = 8;
//   for (int i = 0; i < N; i += TILE_SIZE) {
//     for (int k = 0; k < N; k += TILE_SIZE) {
//       for (int j = 0; j < N; j += TILE_SIZE) {
//         for (int j1 = j; j1 < j + TILE_SIZE; j1++) {
//           for (int i1 = i; i1 < i + TILE_SIZE; i1++) {
//             for (int k1 = k; k1 < k + TILE_SIZE; k1++) {
//               result[i1][j1] += a[i1][k1] * b[k1][j1];
//             }
//           }
//         }
//       }
//     }
//   }
// }

// Compute integer power of a given square matrix
Matrix power(const Matrix &input, const uint32_t k) {
  // Temporary products
  std::unique_ptr<Matrix> productCurrent(new Matrix());
  std::unique_ptr<Matrix> productNext(new Matrix());

  // Temporary elements = a^(2^integer)
  std::unique_ptr<Matrix> elementCurrent(new Matrix());
  std::unique_ptr<Matrix> elementNext(new Matrix());

  // Initial values
  identity(*productCurrent);
  *elementCurrent = input;

  // Use binary representation of k to be O(log(k))
  for (auto i = k; i > 0; i /= 2) {
    if (i % 2 != 0) {
      // Multiply the product by element
      multiply(*productNext, *productCurrent, *elementCurrent);
      std::swap(productNext, productCurrent);

      // Exit early to skip next squaring
      if (i == 1)
        break;
    }

    // Square an element
    multiply(*elementNext, *elementCurrent, *elementCurrent);
    std::swap(elementNext, elementCurrent);
  }

  return std::move(*productCurrent);
}
