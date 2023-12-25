#include "solution.hpp"
#include <algorithm>
#include <cassert>
#include <type_traits>
#include <immintrin.h>

using simd_score_t = std::array<int16_t, sequence_count_v>;
using simd_sequence_t = std::array<simd_score_t, sequence_size_v>;

simd_sequence_t transpond_seq(std::vector<sequence_t> const& seq) {

  simd_sequence_t result;
  for (int i = 0; i < seq.size(); i++ ) {
    for (int j = 0; j < sequence_size_v; j++) {
      result[i][j] = seq[j][i];
    }
  }

  return result;
}

void store_simd_score(simd_score_t& score, __m256i regscore) {
  _mm256_store_si256((__m256i*)score.data(), regscore);
}


__m256i load_simd_score(simd_score_t& score) {
  return _mm256_load_si256((__m256i*)score.data());
}



// The alignment algorithm which computes the alignment of the given sequence
// pairs.
result_t compute_alignment(std::vector<sequence_t> const &sequences1,
                           std::vector<sequence_t> const &sequences2) {
  result_t result{};

 using score_t = int16_t;
 using column_t = std::array<simd_score_t, sequence_size_v + 1>;
 //using column_t = std::array<score_t, sequence_size_v + 1>;

 // sequence_t const &sequence1 = sequences1[sequence_idx];
 // sequence_t const &sequence2 = sequences2[sequence_idx];

 alignas(64)
 auto trn_seq1 = transpond_seq(sequences1);
 alignas(64)
 auto trn_seq2 = transpond_seq(sequences2);

 /*
  * Initialise score values.
  */
 __m256i gap_open = _mm256_set1_epi16(-11);
 __m256i gap_extension = _mm256_set1_epi16(-1);
 __m256i match = _mm256_set1_epi16(6);
 __m256i mismatch = _mm256_set1_epi16(-4);

 /*
  * Setup the matrix.
  * Note we can compute the entire matrix with just one column in memory,
  * since we are only interested in the last value of the last column in the
  * score matrix.
  */
 alignas(64)
 column_t score_column{};
 column_t horizontal_gap_column{};
 __m256i last_vertical_gap{};

 /*
  * Initialise the first column of the matrix.
  */
 // horizontal_gap_column[0] = gap_open;
 store_simd_score(horizontal_gap_column[0] , gap_open);
 last_vertical_gap = gap_open;

 for (size_t i = 1; i < score_column.size(); ++i) {
   store_simd_score(score_column[i] , last_vertical_gap);
   store_simd_score(horizontal_gap_column[i], _mm256_add_epi16( last_vertical_gap , gap_open));
   last_vertical_gap = _mm256_add_epi16(last_vertical_gap, gap_extension);
 }

 /*
  * Compute the main recursion to fill the matrix.
  */
 for (unsigned col = 1; col <= trn_seq2.size(); ++col) {
   // score_t last_diagonal_score = score_column[0]; // Cache last diagonal score to compute this cell.
   // score_column[0] = horizontal_gap_column[0];
   // last_vertical_gap = horizontal_gap_column[0] + gap_open;
   // horizontal_gap_column[0] += gap_extension;

   __m256i last_diagonal_score = load_simd_score(score_column[0]);
   __m256i horizontal_gap_column_simd = load_simd_score(horizontal_gap_column[0]);
   score_column[0] = horizontal_gap_column[0];
   // store_simd_score(score_column[0], horizontal_gap_column_simd);
   last_vertical_gap = _mm256_add_epi16(horizontal_gap_column_simd, gap_open);
   store_simd_score(horizontal_gap_column[0], _mm256_add_epi16(horizontal_gap_column_simd, gap_extension));

   for (unsigned row = 1; row <= trn_seq1.size(); ++row) {
     // Compute next score from diagonal direction with match/mismatch.
     // score_t best_cell_score =
     //    last_diagonal_score +
     //    (sequence1[row - 1] == sequence2[col - 1] ? match : mismatch);

     __m256i seq1row = load_simd_score(trn_seq1[row-1]);
     __m256i seq2col = load_simd_score(trn_seq2[col-1]);
     __m256i eq = _mm256_cmpeq_epi16(seq1row, seq2col);

     __m256i best_cell_score = _mm256_add_epi16(last_diagonal_score, _mm256_and_si256(eq, match));
     best_cell_score = _mm256_add_epi16(best_cell_score, _mm256_andnot_si256(eq, mismatch));

     horizontal_gap_column_simd = load_simd_score(horizontal_gap_column[row]);

     //      // Determine best score from diagonal, vertical, or horizontal
     // // direction.
     // best_cell_score = std::max(best_cell_score, last_vertical_gap);
     // best_cell_score = std::max(best_cell_score, horizontal_gap_column[row]);

     best_cell_score = _mm256_max_epi16(best_cell_score, last_vertical_gap);
     best_cell_score = _mm256_max_epi16(best_cell_score, horizontal_gap_column_simd);

     // // Cache next diagonal value and store optimum in score_column.
     // last_diagonal_score = score_column[row];
     // score_column[row] = best_cell_score;
     last_diagonal_score = load_simd_score(score_column[row]);
     store_simd_score(score_column[row], best_cell_score);

     // // Compute the next values for vertical and horizontal gap.
     // best_cell_score += gap_open;
     // last_vertical_gap += gap_extension;
     // horizontal_gap_column[row] += gap_extension;
      
     best_cell_score = _mm256_add_epi16(best_cell_score, gap_open);
     last_vertical_gap = _mm256_add_epi16(last_vertical_gap, gap_extension);
     horizontal_gap_column_simd = _mm256_add_epi16(horizontal_gap_column_simd, gap_extension);

     // // Store optimum between gap open and gap extension.
     // last_vertical_gap = std::max(last_vertical_gap, best_cell_score);
     // horizontal_gap_column[row] =
     //     std::max(horizontal_gap_column[row], best_cell_score);
     last_vertical_gap = _mm256_max_epi16(last_vertical_gap, best_cell_score);
     store_simd_score(horizontal_gap_column[row], _mm256_max_epi16(horizontal_gap_column_simd, best_cell_score));
   }
 }

 // // Report the best score.
 // result[sequence_idx] = score_column.back();
  

  return score_column.back();
}

// result_t compute_alignment(std::vector<sequence_t> const &sequences1,
//                            std::vector<sequence_t> const &sequences2) {
//   result_t result{};
// 
//   for (size_t sequence_idx = 0; sequence_idx < sequences1.size();
//        ++sequence_idx) {
//     using score_t = int16_t;
//     using column_t = std::array<score_t, sequence_size_v + 1>;
// 
//     sequence_t const &sequence1 = sequences1[sequence_idx];
//     sequence_t const &sequence2 = sequences2[sequence_idx];
// 
//     /*
//      * Initialise score values.
//      */
//     score_t gap_open{-11};
//     score_t gap_extension{-1};
//     score_t match{6};
//     score_t mismatch{-4};
// 
//     /*
//      * Setup the matrix.
//      * Note we can compute the entire matrix with just one column in memory,
//      * since we are only interested in the last value of the last column in the
//      * score matrix.
//      */
//     column_t score_column{};
//     column_t horizontal_gap_column{};
//     score_t last_vertical_gap{};
// 
//     /*
//      * Initialise the first column of the matrix.
//      */
//     horizontal_gap_column[0] = gap_open;
//     last_vertical_gap = gap_open;
// 
//     for (size_t i = 1; i < score_column.size(); ++i) {
//       score_column[i] = last_vertical_gap;
//       horizontal_gap_column[i] = last_vertical_gap + gap_open;
//       last_vertical_gap += gap_extension;
//     }
// 
//     /*
//      * Compute the main recursion to fill the matrix.
//      */
//     for (unsigned col = 1; col <= sequence2.size(); ++col) {
//       score_t last_diagonal_score =
//           score_column[0]; // Cache last diagonal score to compute this cell.
//       score_column[0] = horizontal_gap_column[0];
//       last_vertical_gap = horizontal_gap_column[0] + gap_open;
//       horizontal_gap_column[0] += gap_extension;
// 
//       for (unsigned row = 1; row <= sequence1.size(); ++row) {
//         // Compute next score from diagonal direction with match/mismatch.
//         score_t best_cell_score =
//             last_diagonal_score +
//             (sequence1[row - 1] == sequence2[col - 1] ? match : mismatch);
//         // Determine best score from diagonal, vertical, or horizontal
//         // direction.
//         best_cell_score = std::max(best_cell_score, last_vertical_gap);
//         best_cell_score = std::max(best_cell_score, horizontal_gap_column[row]);
//         // Cache next diagonal value and store optimum in score_column.
//         last_diagonal_score = score_column[row];
//         score_column[row] = best_cell_score;
//         // Compute the next values for vertical and horizontal gap.
//         best_cell_score += gap_open;
//         last_vertical_gap += gap_extension;
//         horizontal_gap_column[row] += gap_extension;
//         // Store optimum between gap open and gap extension.
//         last_vertical_gap = std::max(last_vertical_gap, best_cell_score);
//         horizontal_gap_column[row] =
//             std::max(horizontal_gap_column[row], best_cell_score);
//       }
//     }
// 
//     // Report the best score.
//     result[sequence_idx] = score_column.back();
//   }
// 
//   return result;
// }
