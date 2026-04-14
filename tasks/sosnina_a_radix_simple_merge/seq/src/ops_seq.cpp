#include "sosnina_a_radix_simple_merge/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "sosnina_a_radix_simple_merge/common/include/common.hpp"

namespace sosnina_a_radix_simple_merge {

namespace {

constexpr int kRadixBits = 8;
constexpr int kRadixSize = 1 << kRadixBits;
constexpr int kNumPasses = sizeof(int) / sizeof(uint8_t);
constexpr uint32_t kSignFlip = 0x80000000U;

void RadixSortLSD(std::vector<int> &data, std::vector<int> &buffer) {
  size_t idx = 0;
  for (int elem : data) {
    buffer[idx++] = static_cast<int>(static_cast<uint32_t>(elem) ^ kSignFlip);
  }
  std::swap(data, buffer);

  for (int pass = 0; pass < kNumPasses; ++pass) {
    std::vector<int> count(kRadixSize + 1, 0);

    for (auto elem : data) {
      auto digit = static_cast<uint8_t>((static_cast<uint32_t>(elem) >> (pass * kRadixBits)) & 0xFF);
      ++count[digit + 1];
    }

    for (int i = 1; i <= kRadixSize; ++i) {
      count[i] += count[i - 1];
    }

    for (auto elem : data) {
      auto digit = static_cast<uint8_t>((static_cast<uint32_t>(elem) >> (pass * kRadixBits)) & 0xFF);
      buffer[count[digit]++] = elem;
    }

    std::swap(data, buffer);
  }

  for (int &elem : data) {
    elem = static_cast<int>(static_cast<uint32_t>(elem) ^ kSignFlip);
  }
}

void SimpleMerge(const std::vector<int> &left, const std::vector<int> &right, std::vector<int> &result) {
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

}  // namespace

SosninaATestTaskSEQ::SosninaATestTaskSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool SosninaATestTaskSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool SosninaATestTaskSEQ::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool SosninaATestTaskSEQ::RunImpl() {
  std::vector<int> &data = GetOutput();
  if (data.size() <= 1) {
    return true;
  }

  size_t mid = data.size() / 2;
  std::vector<int> left_part(data.begin(), data.begin() + static_cast<std::ptrdiff_t>(mid));
  std::vector<int> right_part(data.begin() + static_cast<std::ptrdiff_t>(mid), data.end());
  std::vector<int> left_buffer(left_part.size());
  std::vector<int> right_buffer(right_part.size());

  RadixSortLSD(left_part, left_buffer);
  RadixSortLSD(right_part, right_buffer);

  SimpleMerge(left_part, right_part, data);

  return std::ranges::is_sorted(data);
}

bool SosninaATestTaskSEQ::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace sosnina_a_radix_simple_merge
