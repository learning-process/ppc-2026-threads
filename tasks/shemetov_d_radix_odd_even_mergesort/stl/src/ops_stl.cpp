#include "shemetov_d_radix_odd_even_mergesort/stl/include/ops_stl.hpp"

#include <algorithm>
#include <climits>
#include <cstddef>
#include <future>
#include <vector>

#include "shemetov_d_radix_odd_even_mergesort/common/include/common.hpp"
#include "util/include/util.hpp"

namespace shemetov_d_radix_odd_even_mergesort {

ShemetovDRadixOddEvenMergeSortSTL::ShemetovDRadixOddEvenMergeSortSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

void ShemetovDRadixOddEvenMergeSortSTL::RadixSort(std::vector<int> &array, size_t left, size_t right) {
  if (left >= right) {
    return;
  }

  int maximum = *std::max_element(array.begin() + static_cast<int>(left), array.begin() + static_cast<int>(right) + 1);

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

void ShemetovDRadixOddEvenMergeSortSTL::OddEvenMerge(std::vector<int> &array, size_t start_offset, size_t segment) {
  if (segment <= 1) {
    return;
  }

  size_t padding = segment / 2;

  for (size_t index = 0; index < padding; index += 1) {
    if (array[start_offset + index] > array[start_offset + padding + index]) {
      std::swap(array[start_offset + index], array[start_offset + padding + index]);
    }
  }

  for (padding = segment / 4; padding > 0; padding /= 2) {
    size_t step = padding * 2;
    size_t num_indices = ((segment - padding) / step) * padding;

    for (size_t index = 0; index < num_indices; index += 1) {
      size_t i = index % padding;
      size_t start_position = padding + ((index / padding) * step);

      if (array[start_offset + start_position + i] > array[start_offset + start_position + i + padding]) {
        std::swap(array[start_offset + start_position + i], array[start_offset + start_position + i + padding]);
      }
    }
  }
}

bool ShemetovDRadixOddEvenMergeSortSTL::ValidationImpl() {
  const auto &[size, array] = GetInput();
  return size > 0 && static_cast<size_t>(size) == array.size();
}

bool ShemetovDRadixOddEvenMergeSortSTL::PreProcessingImpl() {
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

bool ShemetovDRadixOddEvenMergeSortSTL::RunImpl() {
  if (power_ <= 1) {
    return true;
  }

  size_t threads = ppc::util::GetNumThreads();

  size_t limit = 1;
  while (limit * 2 <= std::min(threads, power_)) {
    limit *= 2;
  }

  size_t chunk_size = power_ / limit;

  std::vector<std::future<void>> sort;
  for (size_t i = 0; i < limit; i += 1) {
    size_t left = i * chunk_size;
    size_t right = left + chunk_size - 1;

    sort.push_back(std::async(std::launch::async, [this, left, right]() {
      ShemetovDRadixOddEvenMergeSortSTL::RadixSort(this->array_, left, right);
    }));
  }
  for (auto &sorted : sort) {
    sorted.get();
  }

  for (size_t segment = chunk_size * 2; segment <= power_; segment *= 2) {
    std::vector<std::future<void>> merge;

    for (size_t i = 0; i < power_; i += segment) {
      merge.push_back(std::async(std::launch::async, [this, i, segment]() {
        ShemetovDRadixOddEvenMergeSortSTL::OddEvenMerge(this->array_, i, segment);
      }));
    }
    for (auto &merged : merge) {
      merged.get();
    }
  }

  return true;
}

bool ShemetovDRadixOddEvenMergeSortSTL::PostProcessingImpl() {
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
