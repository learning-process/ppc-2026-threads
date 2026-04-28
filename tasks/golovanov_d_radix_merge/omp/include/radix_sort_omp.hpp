#pragma once

#include <cstddef>
#include <vector>

class RadixSortOMP {
 public:
  static void Sort(std::vector<double> &arr);

 private:
  static void SortRange(std::vector<double> &arr, std::size_t left, std::size_t right);
  static std::vector<double> Merge(const std::vector<double> &a, const std::vector<double> &b);
};
