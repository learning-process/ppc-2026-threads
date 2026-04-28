#include "popova_e_radix_sort_for_double_with_simple_merge/omp/include/ops_omp.hpp"

#include <omp.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <random>
#include <vector>

#include "popova_e_radix_sort_for_double_with_simple_merge/common/include/common.hpp"

namespace popova_e_radix_sort_for_double_with_simple_merge_threads {

namespace {

uint64_t DoubleToSortable(double value) {
  uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(double));
  bool is_negative = (bits >> 63) == 1;
  if (is_negative) {
    bits = ~bits;
  } else {
    bits ^= (1ULL << 63);
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
  std::vector<uint64_t> buffer(arr.size());

  for (int byte_index = 0; byte_index < bytes_count; byte_index++) {
    int sdvig = byte_index * 8;
    std::array<size_t, base> count = {0};

    for (const auto &val : arr) {
      count.at((val >> sdvig) & 0xFF)++;
    }

    size_t offset = 0;
    for (auto &c : count) {
      size_t tmp = c;
      c = offset;
      offset += tmp;
    }

    for (const auto &val : arr) {
      size_t pos = (val >> sdvig) & 0xFF;
      buffer.at(count.at(pos)) = val;
      count.at(pos)++;
    }
    arr = buffer;
  }
}

std::vector<double> MergeSorted(const std::vector<double> &left, const std::vector<double> &right) {
  std::vector<double> res;
  res.reserve(left.size() + right.size());
  size_t i = 0;
  size_t j = 0;
  while (i < left.size() && j < right.size()) {
    if (left[i] <= right[j]) {
      res.push_back(left[i++]);
    } else {
      res.push_back(right[j++]);
    }
  }
  while (i < left.size()) {
    res.push_back(left[i++]);
  }
  while (j < right.size()) {
    res.push_back(right[j++]);
  }
  return res;
}

double RandomDouble(double min_val, double max_val) {
  static std::random_device rd;
  static std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(min_val, max_val);
  return dis(gen);
}

bool IsSorted(const std::vector<double> &arr) {
  for (size_t i = 1; i < arr.size(); i++) {
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
    std::memcpy(&bits, &value, sizeof(double));
    hash_original ^= bits;
  }

  for (const double &value : result) {
    uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(double));
    hash_result ^= bits;
  }

  return hash_original == hash_result;
}

}  // namespace

PopovaERadixSorForDoubleWithSimpleMergeOMP::PopovaERadixSorForDoubleWithSimpleMergeOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool PopovaERadixSorForDoubleWithSimpleMergeOMP::ValidationImpl() {
  return GetInput() > 0;
}

bool PopovaERadixSorForDoubleWithSimpleMergeOMP::PreProcessingImpl() {
  int size = GetInput();
  array_.resize(size);
  for (int i = 0; i < size; i++) {
    array_[i] = RandomDouble(-100.0, 100.0);
  }
  return true;
}

bool PopovaERadixSorForDoubleWithSimpleMergeOMP::RunImpl() {
  int n_threads = omp_get_max_threads();
  int n = static_cast<int>(array_.size());
  std::vector<std::vector<double>> local_results(n_threads);

  auto &ref_array = array_;
  auto &ref_local_results = local_results;

#pragma omp parallel num_threads(n_threads) default(none) shared(n, n_threads, ref_array, ref_local_results)
  {
    int thread_id = omp_get_thread_num();
    int left_idx = (thread_id * n) / n_threads;
    int right_idx = ((thread_id + 1) * n) / n_threads;

    if (left_idx < right_idx) {
      int local_size = right_idx - left_idx;
      std::vector<uint64_t> local_bits(local_size);

      for (int i = 0; i < local_size; i++) {
        local_bits.at(i) = DoubleToSortable(ref_array.at(left_idx + i));
      }

      RadixSortUInt(local_bits);

      ref_local_results.at(thread_id).resize(local_size);
      for (int i = 0; i < local_size; i++) {
        ref_local_results.at(thread_id).at(i) = SortableToDouble(local_bits.at(i));
      }
    }
  }

  result_.clear();
  if (!local_results.empty()) {
    result_ = local_results.at(0);
    for (int i = 1; i < n_threads; i++) {
      if (!local_results.at(i).empty()) {
        result_ = MergeSorted(result_, local_results.at(i));
      }
    }
  }

  return true;
}

bool PopovaERadixSorForDoubleWithSimpleMergeOMP::PostProcessingImpl() {
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
