#include "shemetov_d_radix_odd_even_mergesort/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <climits>
#include <cstddef>
#include <vector>

#include "oneapi/tbb/info.h"
#include "oneapi/tbb/parallel_for.h"
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
  std::vector<int> position(256);

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

void ShemetovDRadixOddEvenMergeSortTBB::OddEvenMerge(std::vector<int> &array, size_t start_offset, size_t segment) {
  if (segment <= 1) {
    return;
  }

  size_t padding = segment / 2;
  tbb::parallel_for(static_cast<size_t>(0), padding, [&](size_t i) {
    if (array[start_offset + i] > array[start_offset + padding + i]) {
      std::swap(array[start_offset + i], array[start_offset + padding + i]);
    }
  });

  for (padding = segment / 4; padding > 0; padding /= 2) {
    size_t step = padding * 2;
    size_t blocks = (segment - padding) / step;

    tbb::parallel_for(static_cast<size_t>(0), blocks, [&](size_t block) {
      size_t start_position = padding + (block * step);

      for (size_t i = 0; i < padding; i++) {
        if (array[start_offset + start_position + i] > array[start_offset + start_position + i + padding]) {
          std::swap(array[start_offset + start_position + i], array[start_offset + start_position + i + padding]);
        }
      }
    });
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

  if (size_ <= 1) {
    power_ = size_;
    return true;
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

  size_t threads = tbb::info::default_concurrency();
  size_t limit = 1;
  while (limit * 2 <= threads && limit * 2 <= power_) {
    limit *= 2;
  }

  size_t chunk_size = power_ / limit;

  tbb::parallel_for(static_cast<size_t>(0), limit, [&](size_t i) {
    size_t left = i * chunk_size;
    size_t right = left + chunk_size - 1;

    RadixSort(array_, left, right);
  });

  for (size_t segment = chunk_size; segment <= power_; segment *= 2) {
    tbb::parallel_for(static_cast<size_t>(0), power_ / (segment * 2), [&](size_t local) {
      size_t i = local * segment * 2;

      OddEvenMerge(array_, i, segment * 2);
    });
  }

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
