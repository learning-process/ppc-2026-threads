#include "popova_e_radix_sort_for_double_with_simple_merge/seq/include/ops_seq.hpp"

#include <cstdint>
#include <cstring>
#include <random>
#include <vector>

#include "popova_e_radix_sort_for_double_with_simple_merge/common/include/common.hpp"

namespace popova_e_radix_sort_for_double_with_simple_merge_threads {
namespace {
uint64_t DoubleToSortable(double value) {
  uint64_t bits = 0;
  memcpy(&bits, &value, sizeof(double));

  bool is_negative = (bits >> 63) == 1;
  if (is_negative) {
    bits = ~bits;
  } else {
    bits ^= (1ULL << 63);  // сдвиг старшего бита
  }
  return bits;
}

double SortableToDouble(uint64_t bits) {
  bool is_negative = (bits >> 63) == 1;
  if (is_negative) {
    bits ^= (1ULL << 63);
  } else {
    bits = ~bits;
  }

  double value = 0;
  memcpy(&value, &bits, sizeof(double));
  return value;
}

void RadixSortUInt(std::vector<uint64_t> &arr) {
  if (arr.empty()) {
    return;
  }

  const int bytes_count = 8;
  const int base = 256;

  for (int byte_index = 0; byte_index < bytes_count; byte_index++) {
    int sdvig = byte_index * 8;
    std::vector<std::vector<uint64_t>> buckets(base);

    for (const auto &value : arr) {
      const auto bucket_index = static_cast<size_t>((value >> sdvig) & 0xFFULL);
      buckets[bucket_index].push_back(value);
    }

    size_t pos = 0;
    for (const auto &bucket : buckets) {
      for (const auto &value : bucket) {
        arr[pos++] = value;
      }
    }
  }
}

double RandomDouble(double min_val = -100.0, double max_val = 100.0) {
  static std::random_device rd;
  static std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(min_val, max_val);
  return dis(gen);
}

std::vector<double> MergeSorted(const std::vector<double> &left, const std::vector<double> &right) {
  std::vector<double> result;
  size_t i = 0;
  size_t j = 0;
  while (i < left.size() && j < right.size()) {
    if (left[i] <= right[j]) {
      result.push_back(left[i++]);
    } else {
      result.push_back(right[j++]);
    }
  }
  while (i < left.size()) {
    result.push_back(left[i++]);
  }
  while (j < right.size()) {
    result.push_back(right[j++]);
  }

  return result;
}

bool IsSorted(const std::vector<double> &arr) {
  for (size_t i = 1; i < arr.size(); ++i) {
    if (arr[i - 1] > arr[i]) {
      return false;
    }
  }
  return true;
}

bool SameData(const std::vector<double> &original, const std::vector<double> &result) {
  uint64_t hash_original = 0;
  uint64_t hash_result = 0;

  for (const double &value : original) {
    uint64_t bits = 0;
    memcpy(&bits, &value, sizeof(double));
    hash_original ^= bits;
  }

  for (const double &value : result) {
    uint64_t bits = 0;
    memcpy(&bits, &value, sizeof(double));
    hash_result ^= bits;
  }

  return hash_original == hash_result;
}
}  // namespace

PopovaERadixSorForDoubleWithSimpleMergeSEQ::PopovaERadixSorForDoubleWithSimpleMergeSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool PopovaERadixSorForDoubleWithSimpleMergeSEQ::ValidationImpl() {
  return GetInput() > 0;
}

bool PopovaERadixSorForDoubleWithSimpleMergeSEQ::PreProcessingImpl() {
  int size = GetInput();
  array_.resize(size);

  for (int i = 0; i < size; ++i) {
    array_[i] = RandomDouble();
  }

  return true;
}

bool PopovaERadixSorForDoubleWithSimpleMergeSEQ::RunImpl() {
  if (array_.empty()) {
    return false;
  }

  size_t mid = array_.size() / 2;
  std::vector<uint64_t> left_bits;
  std::vector<uint64_t> right_bits;
  left_bits.reserve(mid);
  right_bits.reserve(array_.size() - mid);

  for (size_t i = 0; i < mid; ++i) {
    left_bits.push_back(DoubleToSortable(array_[i]));
  }
  for (size_t i = mid; i < array_.size(); ++i) {
    right_bits.push_back(DoubleToSortable(array_[i]));
  }

  RadixSortUInt(left_bits);
  RadixSortUInt(right_bits);

  std::vector<double> left(left_bits.size());
  std::vector<double> right(right_bits.size());

  for (size_t i = 0; i < left_bits.size(); ++i) {
    left[i] = SortableToDouble(left_bits[i]);
  }
  for (size_t i = 0; i < right_bits.size(); ++i) {
    right[i] = SortableToDouble(right_bits[i]);
  }

  result_ = MergeSorted(left, right);

  return true;
}

bool PopovaERadixSorForDoubleWithSimpleMergeSEQ::PostProcessingImpl() {
  bool sorted = IsSorted(result_);
  bool same = SameData(array_, result_);

  if (sorted && same) {
    GetOutput() = 1;
  } else {
    GetOutput() = 0;
  }

  return true;
}

}  // namespace popova_e_radix_sort_for_double_with_simple_merge_threads
