#include "solution.hpp"
#include <iostream>
#include <immintrin.h>
#include <bit>


// Find the longest line in a file.
// Implementation uses ternary operator with a hope that compiler will
// turn it into a CMOV instruction.
// The code inside the inner loop is equivalent to:
/*
if (s == '\n') {
  longestLine = std::max(curLineLength, longestLine);
  curLineLength = 0;
} else {
  curLineLength++;
}*/

int ntz(uint64_t x) ;

unsigned solution(const std::string &inputContents) {
  unsigned longestLine = 0;
  unsigned curLineLength = 0;

  auto iterations = inputContents.size() / 32;

  uint32_t curLineBegin = 0;

  for (int i = 1; i < iterations; i++) {
    // curLineLength = (s == '\n') ? 0 : curLineLength + 1;
    // longestLine = std::max(curLineLength, longestLine);
    auto inp = _mm256_loadu_si256((__m256i*)(inputContents.data() + i*32));
    auto data = _mm256_and_si256(inp, _mm256_set1_epi8(0xf));
    auto endl_mask = _mm256_cmpeq_epi8(data, _mm256_set1_epi8(10));
    auto mask = (uint32_t)(_mm256_movemask_epi8(endl_mask));
    if (mask == 0) {
      curLineLength += 32;
    } else {
      unsigned pos_count = 0;
      do {
        int pos = _tzcnt_u32(mask);
        curLineLength += pos ;
        if (curLineLength > longestLine) {
          longestLine = curLineLength;
        }
        curLineLength = 0; 
        mask >>= pos;
        pos_count += pos;
      } while(mask > 0);

      curLineLength = 31 - pos_count;

      //curLineLength = 32 - prev_one - 1;

    }
  }

  for (int i = 0; i < inputContents.size()%32; i++) {
    curLineLength = (inputContents[i] == '\n') ? 0 : curLineLength + 1;
    longestLine = std::max(curLineLength, longestLine);
  }
 
  return longestLine;
}

