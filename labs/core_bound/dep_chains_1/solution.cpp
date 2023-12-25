#include "solution.hpp"
#include <array>
#include <iostream>
#include <unordered_set>

unsigned getSumOfDigits(unsigned n) {
  unsigned sum = 0;
  while (n != 0) {
    sum = sum + n % 10;
    n = n / 10;
  }
  return sum;
}

// Task: lookup all the values from l2 in l1.
// For every found value, find the sum of its digits.
// Return the sum of all digits in every found number.
// Both lists have no duplicates and elements placed in *random* order.
// Do NOT sort any of the lists. Do NOT store elements in a hash_map/sets.

// Hint: Traversing a linked list is a long data dependency chain:
//       to get the node N+1 you need to retrieve the node N first.
//       Think how you can execute multiple dependency chains in parallel.
unsigned solution(List *l1, List *l2) {
  unsigned retVal = 0;

  std::unordered_set<unsigned> values;
  while (l1) {
    values.insert(l1->value);
    l1 = l1->next;
  }

  while (l2) {
    if (values.find(l2->value) != values.end()) {
      retVal += getSumOfDigits(l2->value);
    }
    l2 = l2->next;
  }

  return retVal;
}
