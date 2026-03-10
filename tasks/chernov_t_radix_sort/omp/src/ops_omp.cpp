#include "chernov_t_radix_sort/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "chernov_t_radix_sort/common/include/common.hpp"

namespace chernov_t_radix_sort {

ChernovTRadixSortOMP::ChernovTRadixSortOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ChernovTRadixSortOMP::ValidationImpl() {
  return true;
}

bool ChernovTRadixSortOMP::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

constexpr int kBitsPerDigit = 8;
constexpr int kRadix = 1 << kBitsPerDigit;
constexpr uint32_t kSignMask = 0x80000000U;

void ChernovTRadixSortOMP::RadixSortLSD(std::vector<int> &data) {
  if (data.empty()) {
    return;
  }

  const size_t n = data.size();
  std::vector<uint32_t> temp(n);

#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n; ++i) {
    temp[i] = static_cast<uint32_t>(data[i]) ^ kSignMask;
  }

  std::vector<uint32_t> buffer(n);

  for (int byte_index = 0; byte_index < 4; ++byte_index) {
    const int shift = byte_index * kBitsPerDigit;

    int num_threads = 0;
#pragma omp parallel
    {
#pragma omp single
      num_threads = omp_get_num_threads();
    }

    std::vector<std::vector<int>> local_counts(num_threads, std::vector<int>(kRadix, 0));

#pragma omp parallel for schedule(static)
    for (int t = 0; t < num_threads; ++t) {
      size_t start = (n * t) / num_threads;
      size_t end = (n * (t + 1)) / num_threads;
      for (size_t i = start; i < end; ++i) {
        int digit = static_cast<int>((temp[i] >> shift) & 0xFFU);
        ++local_counts[t][digit];
      }
    }

    std::vector<int> count(kRadix, 0);
    for (int t = 0; t < num_threads; ++t) {
      for (int d = 0; d < kRadix; ++d) {
        count[static_cast<size_t>(d)] += local_counts[t][d];
      }
    }

    for (int i = 1; i < kRadix; ++i) {
      count[static_cast<size_t>(i)] += count[static_cast<size_t>(i - 1)];
    }

    std::vector<std::vector<int>> thread_pos(num_threads, std::vector<int>(kRadix, 0));

    for (int t = 0; t < num_threads; ++t) {
      for (int d = 0; d < kRadix; ++d) {
        if (t == 0) {
          thread_pos[t][d] = (d == 0) ? 0 : count[static_cast<size_t>(d - 1)];
        } else {
          thread_pos[t][d] = thread_pos[t - 1][d] + local_counts[t - 1][d];
        }
      }
    }

#pragma omp parallel for schedule(static)
    for (int t = 0; t < num_threads; ++t) {
      size_t start = (n * t) / num_threads;
      size_t end = (n * (t + 1)) / num_threads;
      for (size_t i = start; i < end; ++i) {
        uint32_t val = temp[i];
        int digit = static_cast<int>((val >> shift) & 0xFFU);
        int pos = thread_pos[t][digit]++;
        buffer[static_cast<size_t>(pos)] = val;
      }
    }

    temp.swap(buffer);
  }

#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n; ++i) {
    data[i] = static_cast<int>(temp[i] ^ kSignMask);
  }
}

void ChernovTRadixSortOMP::SimpleMerge(const std::vector<int> &left, const std::vector<int> &right,
                                       std::vector<int> &result) {
  result.resize(left.size() + right.size());

  size_t i = 0, j = 0, k = 0;

  while (i < left.size() && j < right.size()) {
    if (left[i] <= right[j]) {
      result[k++] = left[i++];
    } else {
      result[k++] = right[j++];
    }
  }
  while (i < left.size()) {
    result[k++] = left[i++];
  }
  while (j < right.size()) {
    result[k++] = right[j++];
  }
}

bool ChernovTRadixSortOMP::RunImpl() {
  auto &data = GetOutput();

  if (data.size() <= 1) {
    return true;
  }

  const size_t mid = data.size() / 2;
  std::vector<int> left(data.begin(), data.begin() + static_cast<std::ptrdiff_t>(mid));
  std::vector<int> right(data.begin() + static_cast<std::ptrdiff_t>(mid), data.end());

  RadixSortLSD(left);
  RadixSortLSD(right);

  SimpleMerge(left, right, data);

  return true;
}

bool ChernovTRadixSortOMP::PostProcessingImpl() {
  return std::is_sorted(GetOutput().begin(), GetOutput().end());
}

}  // namespace chernov_t_radix_sort
