#include "popova_e_radix_sort_for_double_with_simple_merge/seq/include/ops_seq.hpp"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <vector>

#include "popova_e_radix_sort_for_double_with_simple_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace popova_e_radix_sort_for_double_with_simple_merge_threads {
namespace {
uint64_t double_to_sortable(double value) {
  uint64_t bits;
  memcpy(&bits, &value, sizeof(double));

  bool is_negative = (bits >> 63) == 1;
  if (is_negative) {
    bits = ~bits;
  } else {
    bits = bits ^ (1ULL << 63);  // сдвиг старшего бита
  }
  return bits;
}

double sortable_to_double(uint64_t bits) {
  bool is_negative = (bits >> 63) == 1;
  if (is_negative) {
    bits = bits ^ (1ULL << 63);
  } else {
    bits = ~bits;
  }

  double value;
  memcpy(&value, &bits, sizeof(double));
  return value;
}

void radix_sort_uint(std::vector<uint64_t> &arr) {
  if (arr.empty()) {
    return;
  }

  const int bytes_count = 8;
  const int base = 256;

  for (int byte_index = 0; byte_index < bytes_count; byte_index++) {
    int sdvig = byte_index * 8;
    std::vector<std::vector<uint64_t>> buckets(base);

    for (size_t i = 0; i < arr.size(); i++) {
      int bucket_index = (arr[i] >> sdvig) & 0xFF;
      buckets[bucket_index].push_back(arr[i]);
    }

    size_t pos = 0;
    for (int b = 0; b < base; b++) {
      for (size_t j = 0; j < buckets[b].size(); j++) {
        arr[pos] = buckets[b][j];
        pos++;
      }
    }
  }
}

double random_double(double min = -100.0, double max = 100.0) {
  int r = rand() % 10000;
  double value = min + (r / 10000.0) * (max - min);
  return value;
}

std::vector<double> merge_sorted(const std::vector<double> &left, const std::vector<double> &right) {
  std::vector<double> result;
  size_t i = 0, j = 0;
  while (i < left.size() && j < right.size()) {
    if (left[i] <= right[j]) {
      result.push_back(left[i]);
      i++;
    } else {
      result.push_back(right[j]);
      j++;
    }
  }
  while (i < left.size()) {
    result.push_back(left[i]);
    i++;
  }
  while (j < right.size()) {
    result.push_back(right[j]);
    j++;
  }
  return result;
}

bool is_sorted(const std::vector<double> &arr) {
  for (size_t i = 1; i < arr.size(); i++) {
    if (arr[i - 1] > arr[i]) {
      return false;
    }
  }
  return true;
}

bool same_data(const std::vector<double> &original, const std::vector<double> &result) {
  uint64_t hash_original = 0;
  uint64_t hash_result = 0;

  for (size_t i = 0; i < original.size(); i++) {
    uint64_t bits;
    memcpy(&bits, &original[i], sizeof(double));
    hash_original ^= bits;
  }
  for (size_t i = 0; i < result.size(); i++) {
    uint64_t bits;
    memcpy(&bits, &result[i], sizeof(double));
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
  if (GetInput() <= 0) {
    return false;
  }
  return true;
}

bool PopovaERadixSorForDoubleWithSimpleMergeSEQ::PreProcessingImpl() {
  int size = GetInput();
  array_.resize(size);

  for (int i = 0; i < size; i++) {
    array_[i] = random_double();
  }

  // std::cout << "[DEBUG] Input array: ";
  // for (size_t i = 0; i < array_.size(); i++) {
  //   std::cout << std::fixed << std::setprecision(2) << array_[i] << " ";
  // }
  // std::cout << std::endl;

  return true;
}

bool PopovaERadixSorForDoubleWithSimpleMergeSEQ::RunImpl() {
  if (array_.empty()) {
    return false;
  }

  size_t mid = array_.size() / 2;
  std::vector<uint64_t> left_bits, right_bits;

  for (size_t i = 0; i < mid; i++) {
    left_bits.push_back(double_to_sortable(array_[i]));
  }
  for (size_t i = mid; i < array_.size(); i++) {
    right_bits.push_back(double_to_sortable(array_[i]));
  }

  radix_sort_uint(left_bits);
  radix_sort_uint(right_bits);

  std::vector<double> left(left_bits.size()), right(right_bits.size());
  for (size_t i = 0; i < left_bits.size(); i++) {
    left[i] = sortable_to_double(left_bits[i]);
  }
  for (size_t i = 0; i < right_bits.size(); i++) {
    right[i] = sortable_to_double(right_bits[i]);
  }

  result_ = merge_sorted(left, right);

  return true;
}

bool PopovaERadixSorForDoubleWithSimpleMergeSEQ::PostProcessingImpl() {
  bool sorted = is_sorted(result_);
  bool same = same_data(array_, result_);

  // std::cout << "[DEBUG] Sorted array: ";
  // for (size_t i = 0; i < result_.size(); i++) {
  //   std::cout << std::fixed << std::setprecision(2) << result_[i] << " ";
  // }
  // std::cout << std::endl;

  if (sorted && same) {
    GetOutput() = 1;
  } else {
    GetOutput() = 0;
  }

  return true;
}

}  // namespace popova_e_radix_sort_for_double_with_simple_merge_threads