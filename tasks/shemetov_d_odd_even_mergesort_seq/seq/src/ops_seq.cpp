#include "shemetov_d_odd_even_mergesort_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <climits>
#include <cstddef>
#include <vector>

#include "shemetov_d_odd_even_mergesort_seq/common/include/common.hpp"

namespace shemetov_d_odd_even_mergesort {

ShemetovDOddEvenMergeSortSEQ::ShemetovDOddEvenMergeSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

void ShemetovDOddEvenMergeSortSEQ::RadixSort(std::vector<int> &array, size_t left, size_t right) {
  if (left >= right) {
    return;
  }

  int maximum = *std::ranges::max_element(array.begin(), array.end());

  size_t segment = right - left + 1;

  this->buffer_.resize(segment);
  for (size_t merge_shift = 0; merge_shift < 32; merge_shift += 8) {
    this->position_.assign(256, 0);

    for (size_t i = left; i <= right; i += 1) {
      int apply_bitmask = (array[i] >> merge_shift) & 0xFF;

      this->position_[apply_bitmask] += 1;
    }

    for (size_t i = 1; i < 256; i += 1) {
      this->position_[i] += this->position_[i - 1];
    }

    for (size_t i = segment; i > 0; i -= 1) {
      size_t current_index = left + i - 1;
      int apply_bitmask = (array[current_index] >> merge_shift) & 0xFF;

      this->buffer_[this->position_[apply_bitmask] -= 1] = array[current_index];
    }

    for (size_t i = 0; i < segment; i += 1) {
      array[left + i] = this->buffer_[i];
    }

    if ((maximum >> merge_shift) < 256) {
      break;
    }
  }
}

void ShemetovDOddEvenMergeSortSEQ::OddEvenMerge(std::vector<int> &array, size_t start_offset, size_t segment) {
  if (segment <= 1) {
    return;
  }

  size_t padding = segment / 2;
  for (size_t i = 0; i < padding; i += 1) {
    if (array[start_offset + i] > array[start_offset + padding + i]) {
      std::swap(array[start_offset + i], array[start_offset + padding + i]);
    }
  }

  for (padding = segment / 4; padding > 0; padding /= 2) {
    size_t step = padding * 2;

    for (size_t start_position = padding; start_position + padding < segment; start_position += step) {
      for (size_t i = 0; i < padding; i += 1) {
        if (array[start_offset + start_position + i] > array[start_offset + start_position + i + padding]) {
          std::swap(array[start_offset + start_position + i], array[start_offset + start_position + i + padding]);
        }
      }
    }
  }
}

bool ShemetovDOddEvenMergeSortSEQ::ValidationImpl() {
  const auto &[size, array] = GetInput();
  return size > 0 && static_cast<size_t>(size) == array.size();
}

bool ShemetovDOddEvenMergeSortSEQ::PreProcessingImpl() {
  const auto &[size, array] = GetInput();

  if (size == 0) {
    return true;
  }

  array_ = array;

  offset_ = *std::ranges::min_element(array_.begin(), array_.end());
  size_ = static_cast<size_t>(size);
  power_ = 1;

  while (power_ < size_) {
    power_ *= 2;
  }

  for (size_t i = 0; i < size_; i += 1) {
    array_[i] -= offset_;
  }

  if (power_ > size_) {
    array_.resize(power_, INT_MAX);
  }

  return true;
}

bool ShemetovDOddEvenMergeSortSEQ::RunImpl() {
  if (power_ <= 1) {
    return true;
  }

  size_t middle = power_ / 2;
  RadixSort(array_, 0, middle - 1);
  RadixSort(array_, middle, power_ - 1);

  OddEvenMerge(array_, 0, power_);

  return true;
}

bool ShemetovDOddEvenMergeSortSEQ::PostProcessingImpl() {
  if (size_ == 0) {
    return true;
  }

  array_.resize(size_);

  for (size_t i = 0; i < size_; i += 1) {
    array_[i] += offset_;
  }

  if (!std::ranges::is_sorted(array_.begin(), array_.end())) {
    return false;
  }

  GetOutput() = array_;
  return true;
}

}  // namespace shemetov_d_odd_even_mergesort
