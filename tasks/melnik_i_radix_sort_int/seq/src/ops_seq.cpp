#include "melnik_i_radix_sort_int/seq/include/ops_seq.hpp"

#include <algorithm>
#include <array>
#include <ranges>
#include <vector>

#include "melnik_i_radix_sort_int/common/include/common.hpp"

namespace melnik_i_radix_sort_int {

namespace {

constexpr int kBitsPerDigit = 8;
constexpr int kBuckets = 1 << kBitsPerDigit;

}  // namespace

MelnikIRadixSortIntSEQ::MelnikIRadixSortIntSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool MelnikIRadixSortIntSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool MelnikIRadixSortIntSEQ::PreProcessingImpl() {
  GetOutput() = GetInput();
  return !GetOutput().empty();
}

bool MelnikIRadixSortIntSEQ::RunImpl() {
  if (GetOutput().empty()) {
    return false;
  }
  RadixSort(GetOutput());
  return !GetOutput().empty();
}

bool MelnikIRadixSortIntSEQ::PostProcessingImpl() {
  return std::ranges::is_sorted(GetOutput());
}

int MelnikIRadixSortIntSEQ::GetMaxValue(const OutType &data) {
  return *std::ranges::max_element(data);
}

void MelnikIRadixSortIntSEQ::CountingSort(OutType &data, int exp, int offset) {
  const auto n = static_cast<int>(data.size());
  if (n == 0) {
    return;
  }

  OutType output(n);
  std::array<int, kBuckets> count{};
  count.fill(0);

  for (int i = 0; i < n; i++) {
    int digit = ((data[i] + offset) / exp) % kBuckets;
    count.at(digit)++;
  }

  for (int i = 1; i < kBuckets; i++) {
    count.at(i) += count.at(i - 1);
  }

  for (int i = n - 1; i >= 0; i--) {
    int digit = ((data[i] + offset) / exp) % kBuckets;
    output[count.at(digit) - 1] = data[i];
    count.at(digit)--;
  }

  data = std::move(output);
}

void MelnikIRadixSortIntSEQ::RadixSort(OutType &data) {
  if (data.empty()) {
    return;
  }

  int max_val = GetMaxValue(data);
  int min_val = *std::ranges::min_element(data);

  if (min_val >= 0) {
    for (int exp = 1; max_val / exp > 0; exp <<= kBitsPerDigit) {
      CountingSort(data, exp, 0);
    }
    return;
  }

  int offset = -min_val;
  for (int exp = 1; (max_val + offset) / exp > 0 || (min_val + offset) / exp > 0; exp <<= kBitsPerDigit) {
    CountingSort(data, exp, offset);
  }
}

}  // namespace melnik_i_radix_sort_int
