#include "shemetov_d_odd_even_mergesort_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <climits>
#include <cstddef>
#include <vector>

#include "shemetov_d_odd_even_mergesort_seq/common/include/common.hpp"

namespace shemetov_d_odd_even_mergesort {

void ShemetovDOddEvenMergeSortSEQ::CompareExchange(int &a, int &b) {
  if (a > b) {
    std::swap(a, b);
  }
}

void ShemetovDOddEvenMergeSortSEQ::PerfectUnshuffle(std::vector<int> &array, size_t left, size_t right) {
  size_t middle = (left + right) / 2;

  for (size_t i = left, j = 0; i < right; i += 2, j += 1) {
    this->buffer_[left + j] = array[i];
    this->buffer_[middle + j + 1] = array[i + 1];
  }

  for (size_t i = left; i <= right; i += 1) {
    array[i] = this->buffer_[i];
  }
}

void ShemetovDOddEvenMergeSortSEQ::PerfectShuffle(std::vector<int> &array, size_t left, size_t right) {
  size_t middle = (left + right) / 2;

  for (size_t i = left, j = 0; i < right; i += 2, j += 1) {
    this->buffer_[i] = array[left + j];
    this->buffer_[i + 1] = array[middle + j + 1];
  }

  for (size_t i = left; i <= right; i += 1) {
    array[i] = this->buffer_[i];
  }
}

void ShemetovDOddEvenMergeSortSEQ::Merge(std::vector<int> &array, size_t left, size_t right) {
  size_t segment = right - left + 1;

  if (segment < 2) {
    return;
  }

  if (segment == 2) {
    CompareExchange(array[left], array[right]);
    return;
  }

  for (size_t merge_step = segment; merge_step > 2; merge_step /= 2) {
    for (size_t start = left; start <= right; start += merge_step) {
      PerfectUnshuffle(array, start, start + merge_step - 1);
    }
  }

  for (size_t i = left; i <= right; i += 2) {
    CompareExchange(array[i], array[i + 1]);
  }

  for (size_t merge_step = 4; merge_step <= segment; merge_step *= 2) {
    for (size_t start = left; start <= right; start += merge_step) {
      PerfectShuffle(array, start, start + merge_step - 1);

      for (size_t i = start + 1; i < start + merge_step - 1; i += 2) {
        CompareExchange(array[i], array[i + 1]);
      }
    }
  }
}

void ShemetovDOddEvenMergeSortSEQ::OddEvenMergesort(std::vector<int> &array, size_t left, size_t right) {
  size_t segment = right - left + 1;

  for (size_t size = 2; size <= segment; size *= 2) {
    for (size_t start = left; start + size - 1 <= right; start += size) {
      Merge(array, start, start + size - 1);
    }
  }
}

void ShemetovDOddEvenMergeSortSEQ::SortWithPadding(std::vector<int> &array, size_t size, size_t power) {
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

  buffer_.assign(power_, 0);

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
