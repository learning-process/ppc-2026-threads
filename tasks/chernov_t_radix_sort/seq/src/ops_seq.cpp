#include "chernov_t_radix_sort/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "chernov_t_radix_sort/common/include/common.hpp"

namespace chernov_t_radix_sort {

ChernovTRadixSortSEQ::ChernovTRadixSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ChernovTRadixSortSEQ::ValidationImpl() {
  return true;
}

bool ChernovTRadixSortSEQ::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

constexpr int kBitsPerDigit = 8;
constexpr int kRadix = 1 << kBitsPerDigit;
constexpr uint32_t kSignMask = 0x80000000U;

void ChernovTRadixSortSEQ::RadixSortLSD(std::vector<int> &data) {
  if (data.empty()) {
    return;
  }

  std::vector<uint32_t> temp(data.size());
  for (size_t i = 0; i < data.size(); ++i) {
    temp[i] = static_cast<uint32_t>(data[i]) ^ kSignMask;
  }

  std::vector<uint32_t> buffer(temp.size());

  for (int byte_index = 0; byte_index < 4; ++byte_index) {
    std::vector<int> count(kRadix, 0);

    for (uint32_t val : temp) {
      int digit = static_cast<int>((val >> (byte_index * kBitsPerDigit)) & 0xFFU);
      ++count[static_cast<size_t>(digit)];
    }

    for (int i = 1; i < kRadix; ++i) {
      count[static_cast<size_t>(i)] += count[static_cast<size_t>(i - 1)];
    }

    for (int i = static_cast<int>(temp.size()) - 1; i >= 0; --i) {
      uint32_t val = temp[i];
      int digit = static_cast<int>((val >> (byte_index * kBitsPerDigit)) & 0xFFU);
      int pos = --count[static_cast<size_t>(digit)];
      buffer[static_cast<size_t>(pos)] = val;
    }

    temp.swap(buffer);
  }

  for (size_t i = 0; i < data.size(); ++i) {
    data[i] = static_cast<int>(temp[i] ^ kSignMask);
  }
}

void ChernovTRadixSortSEQ::SimpleMerge(const std::vector<int> &left, const std::vector<int> &right,
                                       std::vector<int> &result) {
  result.resize(left.size() + right.size());

  size_t idx_left = 0;
  size_t idx_right = 0;
  size_t idx_result = 0;

  while (idx_left < left.size() && idx_right < right.size()) {
    if (left[idx_left] <= right[idx_right]) {
      result[idx_result++] = left[idx_left++];
    } else {
      result[idx_result++] = right[idx_right++];
    }
  }

  while (idx_left < left.size()) {
    result[idx_result++] = left[idx_left++];
  }
  while (idx_right < right.size()) {
    result[idx_result++] = right[idx_right++];
  }
}

bool ChernovTRadixSortSEQ::RunImpl() {
  auto &data = GetOutput();

  if (data.size() <= 1) {
    return true;
  }

  size_t middle_index = data.size() / 2;
  std::vector<int> left_part(data.begin(), data.begin() + static_cast<std::ptrdiff_t>(middle_index));
  std::vector<int> right_part(data.begin() + static_cast<std::ptrdiff_t>(middle_index), data.end());

  RadixSortLSD(left_part);
  RadixSortLSD(right_part);

  SimpleMerge(left_part, right_part, data);

  return true;
}

bool ChernovTRadixSortSEQ::PostProcessingImpl() {
  return std::is_sorted(GetOutput().begin(), GetOutput().end());
}

}  // namespace chernov_t_radix_sort
