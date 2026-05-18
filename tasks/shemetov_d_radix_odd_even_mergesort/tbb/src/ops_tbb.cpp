#include "shemetov_d_radix_odd_even_mergesort/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>
#include <tbb/tbb.h>

#include <algorithm>
#include <climits>
#include <cstddef>
#include <vector>

#include "shemetov_d_radix_odd_even_mergesort/common/include/common.hpp"

namespace shemetov_d_radix_odd_even_mergesort {

ShemetovDRadixOddEvenMergeSortTBB::ShemetovDRadixOddEvenMergeSortTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

void ShemetovDRadixOddEvenMergeSortTBB::RadixSort(std::vector<int> &array, size_t left, size_t right) {
  if (left >= right) {
    return;
  }

  int maximum =
      *std::ranges::max_element(array.begin() + static_cast<int>(left), array.begin() + static_cast<int>(right) + 1);

  size_t segment = right - left + 1;

  std::vector<int> buffer(segment);
  std::vector<int> position(256, 0);

  for (size_t merge_shift = 0; merge_shift < 32; merge_shift += 8) {
    position.assign(256, 0);

    for (size_t i = left; i <= right; i += 1) {
      int apply_bitmask = (array[i] >> merge_shift) & 0xFF;

      position[apply_bitmask] += 1;
    }

    for (size_t i = 1; i < 256; i += 1) {
      position[i] += position[i - 1];
    }

    for (size_t i = segment; i > 0; i -= 1) {
      size_t current_index = left + i - 1;
      int apply_bitmask = (array[current_index] >> merge_shift) & 0xFF;

      buffer[position[apply_bitmask] -= 1] = array[current_index];
    }

    for (size_t i = 0; i < segment; i += 1) {
      array[left + i] = buffer[i];
    }

    if ((maximum >> merge_shift) < 256) {
      break;
    }
  }
}

void ShemetovDRadixOddEvenMergeSortTBB::ApplyFirstPass(std::vector<int> &array, size_t start_offset, size_t padding) {
  tbb::parallel_for(tbb::blocked_range<size_t>(0, padding), [&](const tbb::blocked_range<size_t> &range) {
    for (size_t i = range.begin(); i != range.end(); ++i) {
      if (array[start_offset + i] > array[start_offset + padding + i]) {
        std::swap(array[start_offset + i], array[start_offset + padding + i]);
      }
    }
  });
}

void ShemetovDRadixOddEvenMergeSortTBB::ApplyMainPass(std::vector<int> &array, size_t start_offset, size_t segment,
                                                      size_t padding) {
  size_t step = padding * 2;
  size_t num_blocks = (segment - padding) / step;

  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_blocks), [&](const tbb::blocked_range<size_t> &range) {
    for (size_t block = range.begin(); block != range.end(); ++block) {
      size_t start_position = padding + (block * step);

      for (size_t i = 0; i < padding; i += 1) {
        if (array[start_offset + start_position + i] > array[start_offset + start_position + i + padding]) {
          std::swap(array[start_offset + start_position + i], array[start_offset + start_position + i + padding]);
        }
      }
    }
  });
}

void ShemetovDRadixOddEvenMergeSortTBB::OddEvenMerge(std::vector<int> &array, size_t start_offset, size_t segment) {
  if (segment <= 1) {
    return;
  }

  size_t padding = segment / 2;
  ApplyFirstPass(array, start_offset, padding);

  for (padding = segment / 4; padding > 0; padding /= 2) {
    ApplyMainPass(array, start_offset, segment, padding);
  }
}

bool ShemetovDRadixOddEvenMergeSortTBB::ValidationImpl() {
  const auto &[size, array] = GetInput();
  return size > 0 && static_cast<size_t>(size) == array.size();
}

bool ShemetovDRadixOddEvenMergeSortTBB::PreProcessingImpl() {
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

bool ShemetovDRadixOddEvenMergeSortTBB::RunImpl() {
  if (power_ <= 1) {
    return true;
  }

  size_t middle = power_ / 2;

  tbb::parallel_invoke([&]() { RadixSort(array_, 0, middle - 1); }, [&]() { RadixSort(array_, middle, power_ - 1); });

  OddEvenMerge(array_, 0, power_);

  return true;
}

bool ShemetovDRadixOddEvenMergeSortTBB::PostProcessingImpl() {
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

}  // namespace shemetov_d_radix_odd_even_mergesort
