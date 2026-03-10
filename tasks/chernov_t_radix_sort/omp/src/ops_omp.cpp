#include "chernov_t_radix_sort/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
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

  for (size_t i = 0; i < n; ++i) {
    temp[i] = static_cast<uint32_t>(data[i]) ^ kSignMask;
  }

  std::vector<uint32_t> buffer(n);

  for (int byte_index = 0; byte_index < 4; ++byte_index) {
    const int shift = byte_index * kBitsPerDigit;

    std::vector<int> count(kRadix, 0);
    for (size_t i = 0; i < n; ++i) {
      int digit = static_cast<int>((temp[i] >> shift) & 0xFFU);
      ++count[static_cast<size_t>(digit)];
    }

    for (int i = 1; i < kRadix; ++i) {
      count[static_cast<size_t>(i)] += count[static_cast<size_t>(i - 1)];
    }

    for (size_t i = n; i-- > 0;) {
      uint32_t val = temp[i];
      int digit = static_cast<int>((val >> shift) & 0xFFU);
      buffer[static_cast<size_t>(--count[static_cast<size_t>(digit)])] = val;
    }

    temp.swap(buffer);
  }

  for (size_t i = 0; i < n; ++i) {
    data[i] = static_cast<int>(temp[i] ^ kSignMask);
  }
}

void ChernovTRadixSortOMP::SimpleMerge(const std::vector<int> &left, const std::vector<int> &right,
                                       std::vector<int> &result) {
  result.resize(left.size() + right.size());

  size_t i = 0;
  size_t j = 0;
  size_t k = 0;

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

#pragma omp parallel sections default(none) shared(left, right, data)
  {
#pragma omp section
    {
      RadixSortLSD(left);
    }

#pragma omp section
    {
      RadixSortLSD(right);
    }
  }

  SimpleMerge(left, right, data);

  return true;
}

bool ChernovTRadixSortOMP::PostProcessingImpl() {
  return std::is_sorted(GetOutput().begin(), GetOutput().end());
}

}  // namespace chernov_t_radix_sort
