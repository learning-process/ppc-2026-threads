#include "shemetov_d_odd_even_mergesort_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <climits>
#include <cstddef>
#include <vector>

#include "shemetov_d_odd_even_mergesort_seq/common/include/common.hpp"

namespace shemetov_d_odd_even_mergesort {

namespace {

void CompareExchange(int &a, int &b) {
  if (a > b) {
    std::swap(a, b);
  }
}

void PerfectUnshuffle(std::vector<int> &array, size_t left, size_t right) {
  size_t middle = (left + right) / 2;
  std::vector<int> temp_array(array.size());

  for (size_t i = left, j = 0; i < right; i += 2, j += 1) {
    temp_array[left + j] = array[i];
    temp_array[middle + j + 1] = array[i + 1];
  }

  for (size_t i = left; i <= right; i += 1) {
    array[i] = temp_array[i];
  }
}

void PerfectShuffle(std::vector<int> &array, size_t left, size_t right) {
  size_t middle = (left + right) / 2;
  std::vector<int> temp_array(array.size());

  for (size_t i = left, j = 0; i < right; i += 2, j += 1) {
    temp_array[i] = array[left + j];
    temp_array[i + 1] = array[middle + j + 1];
  }

  for (size_t i = left; i <= right; i += 1) {
    array[i] = temp_array[i];
  }
}

void Merge(std::vector<int> &array, size_t left, size_t right) {
  if (right == left + 1) {
    CompareExchange(array[left], array[right]);
    return;
  }

  if (right < left + 2) {
    return;
  }

  size_t middle = (left + right) / 2;

  PerfectUnshuffle(array, left, right);

  Merge(array, left, middle);
  Merge(array, middle + 1, right);

  PerfectShuffle(array, left, right);

  for (size_t i = left + 1; i < right; i += 2) {
    CompareExchange(array[i], array[i + 1]);
  }
}

void OddEvenMergesort(std::vector<int> &array, size_t left, size_t right) {
  if (right <= left) {
    return;
  }

  size_t middle = left + ((right - left) / 2);

  OddEvenMergesort(array, left, middle);
  OddEvenMergesort(array, middle + 1, right);

  Merge(array, left, right);
}

void SortWithPadding(std::vector<int> &array, size_t size, size_t power) {
  if (power == size) {
    OddEvenMergesort(array, 0, power - 1);
    return;
  }

  std::vector<int> temp_array(power, INT_MAX);

  for (size_t i = 0; i < size; i += 1) {
    temp_array[i] = array[i];
  }

  OddEvenMergesort(temp_array, 0, power - 1);

  for (size_t i = 0; i < size; i += 1) {
    array[i] = temp_array[i];
  }
}

}  // namespace

ShemetovDOddEvenMergeSortSEQ::ShemetovDOddEvenMergeSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ShemetovDOddEvenMergeSortSEQ::ValidationImpl() {
  const auto &[size, array] = GetInput();

  return size > 0 && static_cast<size_t>(size) == array.size();
}

bool ShemetovDOddEvenMergeSortSEQ::PreProcessingImpl() {
  const auto &[size, array] = GetInput();

  array_ = array;

  size_ = static_cast<size_t>(size);
  power_ = 1;

  while (power_ < size_) {
    power_ *= 2;
  }

  return true;
}

bool ShemetovDOddEvenMergeSortSEQ::RunImpl() {
  SortWithPadding(array_, size_, power_);

  return true;
}

bool ShemetovDOddEvenMergeSortSEQ::PostProcessingImpl() {
  if (!std::ranges::is_sorted(array_.begin(), array_.end())) {
    return false;
  }

  GetOutput() = array_;

  return true;
}

}  // namespace shemetov_d_odd_even_mergesort
