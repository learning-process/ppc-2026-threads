#include "votincev_d_radixmerge_sort/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <future>
#include <utility>
#include <vector>

#include "votincev_d_radixmerge_sort/common/include/common.hpp"

namespace votincev_d_radixmerge_sort {

VotincevDRadixMergeSortSTL::VotincevDRadixMergeSortSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool VotincevDRadixMergeSortSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool VotincevDRadixMergeSortSTL::PreProcessingImpl() {
  return true;
}

void VotincevDRadixMergeSortSTL::LocalRadixSort(uint32_t *begin, uint32_t *end) {
  auto n = static_cast<int32_t>(end - begin);
  if (n <= 1) {
    return;
  }

  uint32_t max_val = *std::max_element(begin, end);

  std::vector<uint32_t> buffer(static_cast<size_t>(n));
  uint32_t *src = begin;
  uint32_t *dst = buffer.data();

  for (int64_t exp = 1; (static_cast<int64_t>(max_val) / exp) > 0; exp *= 10) {
    std::array<int32_t, 10> count{};

    for (int32_t i = 0; i < n; ++i) {
      count.at(static_cast<size_t>((src[i] / exp) % 10))++;
    }
    for (int32_t i = 1; i < 10; ++i) {
      count.at(static_cast<size_t>(i)) += count.at(static_cast<size_t>(i - 1));
    }
    for (int32_t i = n - 1; i >= 0; --i) {
      auto digit = static_cast<size_t>((src[i] / exp) % 10);
      size_t target_idx = static_cast<size_t>(count.at(digit)) - 1;
      dst[target_idx] = src[i];
      count.at(digit)--;
    }
    std::swap(src, dst);
  }

  if (src != begin) {
    std::copy(src, src + n, begin);
  }
}

void VotincevDRadixMergeSortSTL::Merge(const uint32_t *src, uint32_t *dst, int32_t left, int32_t mid, int32_t right) {
  int32_t i = left;
  int32_t j = mid;
  int32_t k = left;
  while (i < mid && j < right) {
    dst[k++] = (src[i] <= src[j]) ? src[i++] : src[j++];
  }
  while (i < mid) {
    dst[k++] = src[i++];
  }
  while (j < right) {
    dst[k++] = src[j++];
  }
}

void VotincevDRadixMergeSortSTL::ParallelRadixMergeSort(uint32_t *data, int32_t n, uint32_t *temp) {
  const int32_t grain_size = 4096;
  if (n <= grain_size) {
    LocalRadixSort(data, data + n);
    return;
  }

  std::vector<std::future<void>> futures;

  int32_t num_blocks = (n + grain_size - 1) / grain_size;
  futures.reserve(static_cast<size_t>(num_blocks));

  for (int32_t i = 0; i < num_blocks; ++i) {
    int32_t left = i * grain_size;
    int32_t right = std::min(left + grain_size, n);

    futures.push_back(std::async(std::launch::async, [data, left, right] {
      VotincevDRadixMergeSortSTL::LocalRadixSort(data + left, data + right);
    }));
  }
  for (auto &f : futures) {
    f.get();
  }
  futures.clear();

  uint32_t *src = data;
  uint32_t *dst = temp;

  for (int32_t width = grain_size; width < n; width *= 2) {
    int32_t num_merges = (n + (2 * width) - 1) / (2 * width);

    for (int32_t i = 0; i < num_merges; ++i) {
      int32_t left = i * 2 * width;
      int32_t mid = std::min(left + width, n);
      int32_t right = std::min(left + 2 * width, n);

      if (mid < right) {
        futures.push_back(std::async(std::launch::async, [src, dst, left, mid, right] {
          VotincevDRadixMergeSortSTL::Merge(src, dst, left, mid, right);
        }));
      } else if (left < n) {
        std::copy(src + left, src + mid, dst + left);
      }
    }
    for (auto &f : futures) {
      f.get();
    }
    futures.clear();

    std::swap(src, dst);
  }

  if (src != data) {
    std::copy(src, src + n, data);
  }
}

bool VotincevDRadixMergeSortSTL::RunImpl() {
  const auto &input = GetInput();
  auto n = static_cast<int32_t>(input.size());
  if (n == 0) {
    return true;
  }

  int32_t min_val = *std::min_element(input.begin(), input.end());

  std::vector<uint32_t> working_array(static_cast<size_t>(n));
  for (int32_t i = 0; i < n; ++i) {
    working_array[static_cast<size_t>(i)] =
        static_cast<uint32_t>(input[static_cast<size_t>(i)]) - static_cast<uint32_t>(min_val);
  }

  std::vector<uint32_t> temp_buffer(static_cast<size_t>(n));
  ParallelRadixMergeSort(working_array.data(), n, temp_buffer.data());

  std::vector<int32_t> result(static_cast<size_t>(n));
  for (int32_t i = 0; i < n; ++i) {
    result[static_cast<size_t>(i)] =
        static_cast<int32_t>(working_array[static_cast<size_t>(i)] + static_cast<uint32_t>(min_val));
  }

  GetOutput() = std::move(result);
  return true;
}

bool VotincevDRadixMergeSortSTL::PostProcessingImpl() {
  return true;
}

}  // namespace votincev_d_radixmerge_sort
