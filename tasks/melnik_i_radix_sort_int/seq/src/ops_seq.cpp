#include "melnik_i_radix_sort_int/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace melnik_i_radix_sort_int {

namespace {

constexpr int kBitsPerDigit = 8;
constexpr int kBuckets = 256;

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
  return std::is_sorted(GetOutput().begin(), GetOutput().end());
}

int MelnikIRadixSortIntSEQ::GetMaxValue(const OutType &data) {
  return *std::max_element(data.begin(), data.end());
}

void MelnikIRadixSortIntSEQ::CountingSort(OutType &data, int exp) {
  const auto n = static_cast<int>(data.size());
  if (n == 0) {
    return;
  }

  OutType output(n);
  int count[kBuckets] = {0};

  for (int i = 0; i < n; i++) {
    int digit = (data[i] / exp) % kBuckets;
    count[digit]++;
  }

  for (int i = 1; i < kBuckets; i++) {
    count[i] += count[i - 1];
  }

  for (int i = n - 1; i >= 0; i--) {
    int digit = (data[i] / exp) % kBuckets;
    output[count[digit] - 1] = data[i];
    count[digit]--;
  }

  data = std::move(output);
}

void MelnikIRadixSortIntSEQ::RadixSort(OutType &data) {
  if (data.empty()) {
    return;
  }

  int max_val = GetMaxValue(data);
  int min_val = *std::min_element(data.begin(), data.end());

  if (min_val >= 0) {
    for (int exp = 1; max_val / exp > 0; exp <<= kBitsPerDigit) {
      CountingSort(data, exp);
    }
    return;
  }

  if (max_val <= 0) {
    for (auto &val : data) {
      val = -val;
    }
    for (int exp = 1; (-min_val) / exp > 0; exp <<= kBitsPerDigit) {
      CountingSort(data, exp);
    }
    for (auto &val : data) {
      val = -val;
    }
    return;
  }

  auto pivot = std::stable_partition(data.begin(), data.end(), [](int x) { return x < 0; });

  OutType negatives(data.begin(), pivot);
  OutType positives(pivot, data.end());

  if (!negatives.empty()) {
    for (auto &val : negatives) {
      val = -val;
    }
    int neg_max = *std::max_element(negatives.begin(), negatives.end());
    for (int exp = 1; neg_max / exp > 0; exp <<= kBitsPerDigit) {
      CountingSort(negatives, exp);
    }
    for (auto &val : negatives) {
      val = -val;
    }
  }

  if (!positives.empty()) {
    int pos_max = *std::max_element(positives.begin(), positives.end());
    for (int exp = 1; pos_max / exp > 0; exp <<= kBitsPerDigit) {
      CountingSort(positives, exp);
    }
  }

  data.clear();
  data.insert(data.end(), negatives.begin(), negatives.end());
  data.insert(data.end(), positives.begin(), positives.end());
}

}  // namespace melnik_i_radix_sort_int
