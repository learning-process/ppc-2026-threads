#include "shemetov_d_radix_odd_even_mergesort/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <climits>
#include <cstddef>
#include <vector>

#include "shemetov_d_radix_odd_even_mergesort/common/include/common.hpp"

namespace shemetov_d_radix_odd_even_mergesort {

ShemetovDRadixOddEvenMergeSortOMP::ShemetovDRadixOddEvenMergeSortOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

void ShemetovDRadixOddEvenMergeSortOMP::RadixSort(std::vector<int> &array, size_t left, size_t right,
                                                  std::vector<int> &buffer, std::vector<int> &position) {
  if (left >= right) {
    return;
  }

  int maximum = *std::ranges::max_element(array.begin(), array.end());

  size_t segment = right - left + 1;

  buffer.resize(segment);
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

void ShemetovDRadixOddEvenMergeSortOMP::OddEvenMerge(std::vector<int> &array, size_t start_offset, size_t segment) {
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

    for (size_t start_position = padding; start_position < segment - padding; start_position += step) {
      for (size_t i = 0; i < padding; i += 1) {
        if (array[start_offset + start_position + i] > array[start_offset + start_position + i + padding]) {
          std::swap(array[start_offset + start_position + i], array[start_offset + start_position + i + padding]);
        }
      }
    }
  }
}

bool ShemetovDRadixOddEvenMergeSortOMP::ValidationImpl() {
  const auto &[size, array] = GetInput();
  return size > 0 && static_cast<size_t>(size) == array.size();
}

bool ShemetovDRadixOddEvenMergeSortOMP::PreProcessingImpl() {
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

bool ShemetovDRadixOddEvenMergeSortOMP::RunImpl() {
  if (power_ <= 1) {
    return true;
  }

  size_t threads = omp_get_max_threads();

  size_t limit = 1;
  while (limit * 2 <= std::min(threads, power_)) {
    limit *= 2;
  }

  size_t chunk_size = power_ / limit;

  std::vector<int> &ref_array = array_;
  size_t ref_power = power_;

#pragma omp parallel num_threads(limit) default(none) shared(ref_array, ref_power, chunk_size, limit)
  {
    size_t thread_num = omp_get_thread_num();
    size_t left = thread_num * chunk_size;
    size_t right = left + chunk_size - 1;

    std::vector<int> buffer;
    std::vector<int> position;

    RadixSort(ref_array, left, right, buffer, position);
#pragma omp barrier
    for (size_t segment = chunk_size * 2; segment <= ref_power; segment *= 2) {
#pragma omp for
      for (size_t i = 0; i < ref_power; i += segment) {
        OddEvenMerge(ref_array, i, segment);
      }
#pragma omp barrier
    }
  }

  return true;
}

bool ShemetovDRadixOddEvenMergeSortOMP::PostProcessingImpl() {
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
