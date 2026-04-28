#pragma once

#include <omp.h>

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

namespace mityaeva_radix {
class SorterOmp {
 private:
  static uint64_t DoubleToSortable(uint64_t x);

  static uint64_t SortableToDouble(uint64_t x);

  static void CountingPass(std::vector<uint64_t> *current, std::vector<uint64_t> *next, int shift, int radix,
                           int num_threads, size_t data_size);

 public:
  static void Sort(std::vector<double> &data);
};
}  // namespace mityaeva_radix
