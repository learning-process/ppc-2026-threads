#include "popova_e_radix_sort_for_double_with_simple_merge/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <random>
#include <vector>

#include "popova_e_radix_sort_for_double_with_simple_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace popova_e_radix_sort_for_double_with_simple_merge_threads {

namespace {

uint64_t DoubleToSortable(double value) {
  uint64_t bits = 0;
  memcpy(&bits, &value, sizeof(double));
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
    size_t count[base] = {0};
    for (size_t i = 0; i < arr.size(); i++) {
      count[(arr[i] >> sdvig) & 0xFF]++;
    }
    size_t offset = 0;
    for (int i = 0; i < base; i++) {
      size_t tmp = count[i];
      count[i] = offset;
      offset += tmp;
    }
    for (size_t i = 0; i < arr.size(); i++) {
      size_t pos = (arr[i] >> sdvig) & 0xFF;
      buffer[count[pos]] = arr[i];
      count[pos]++;
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

  // std::cout << "\n--- [STEP 1] ORIGINAL ARRAY ---" << std::endl;
  // std::cout << "Size: " << n << "\nData: ";
  // for (int i = 0; i < std::min(n, 10); ++i) std::cout << array_[i] << " ";
  // if (n > 10) std::cout << "...";
  // std::cout << std::endl;

  // std::cout << "\n--- [STEP 2] THREAD DISTRIBUTION ---" << std::endl;

#pragma omp parallel num_threads(n_threads)
  {
    int thread_id = omp_get_thread_num();
    int left_idx = (thread_id * n) / n_threads;
    int right_idx = ((thread_id + 1) * n) / n_threads;

    if (left_idx < right_idx) {
      int local_size = right_idx - left_idx;

      // #pragma omp critical
      // {
      //     std::cout << "[Thread " << thread_id << "] Range [" << left_idx << ", " << right_idx
      //               << "), Data: ";
      //     for (int i = 0; i < std::min(local_size, 5); ++i)
      //         std::cout << array_[left_idx + i] << " ";
      //     if (local_size > 5) std::cout << "...";
      //     std::cout << std::endl;
      // }

      std::vector<uint64_t> local_bits(local_size);
      for (int i = 0; i < local_size; i++) {
        local_bits[i] = DoubleToSortable(array_[left_idx + i]);
      }

      RadixSortUInt(local_bits);

      local_results[thread_id].resize(local_size);
      for (int i = 0; i < local_size; i++) {
        local_results[thread_id][i] = SortableToDouble(local_bits[i]);
      }

      // #pragma omp critical
      // {
      //     std::cout << "[Thread " << thread_id << "] SORTED: ";
      //     for (int i = 0; i < std::min(local_size, 5); ++i)
      //         std::cout << local_results[thread_id][i] << " ";
      //     if (local_size > 5) std::cout << "...";
      //     std::cout << std::endl;
      // }
    }
  }

  // std::cout << "\n--- [STEP 3] MERGING ---" << std::endl;
    result_.clear();
    result_ = local_results[0];

  for (int i = 1; i < n_threads; i++) {
    if (!local_results[i].empty()) {
      // std::cout << "[MASTER] Merging with Thread " << i << "..." << std::endl;
      result_ = MergeSorted(result_, local_results[i]);
    }
  }

  // std::cout << "\n--- [STEP 4] FINAL RESULT ---" << std::endl;
  // for (int i = 0; i < std::min(n, 10); ++i) std::cout << result_[i] << " ";
  // if (n > 10) std::cout << "...";
  // std::cout << "\n-----------------------------\n" << std::endl;

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
