#pragma once

#include <vector>

namespace mityaeva_radix {
class SorterSeq {
 public:
  static void CountingSortAsc(std::vector<double> &source, std::vector<double> &destination, int byte);
  static void CountingSortDesc(std::vector<double> &source, std::vector<double> &destination, int byte);
  static void LSDSortDouble(std::vector<double> &inp);
};
}  // namespace mityaeva_radix
