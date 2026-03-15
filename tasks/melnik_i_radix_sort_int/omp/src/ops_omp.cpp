#include "melnik_i_radix_sort_int/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <vector>

namespace melnik_i_radix_sort_int {

namespace {

constexpr int kBitsPerDigit = 8;
constexpr int kBuckets = 1 << kBitsPerDigit;

}  // namespace

MelnikIRadixSortIntOMP::MelnikIRadixSortIntOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool MelnikIRadixSortIntOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool MelnikIRadixSortIntOMP::PreProcessingImpl() {
  GetOutput() = GetInput();
  return !GetOutput().empty();
}

bool MelnikIRadixSortIntOMP::RunImpl() {
  if (GetOutput().empty()) {
    return false;
  }
  RadixSort(GetOutput());
  return !GetOutput().empty();
}

bool MelnikIRadixSortIntOMP::PostProcessingImpl() {
  return std::is_sorted(GetOutput().begin(), GetOutput().end());
}

int MelnikIRadixSortIntOMP::GetMaxValue(const OutType &data) {
  return *std::max_element(data.begin(), data.end());
}

void MelnikIRadixSortIntOMP::ParallelCountingSort(OutType &data, int exp, int offset) {
  const auto n = static_cast<int>(data.size());
  if (n == 0) {
    return;
  }

  const int num_threads = omp_get_max_threads();
  const int chunk_size = (n + num_threads - 1) / num_threads;

  std::vector<int> local_counts(num_threads * kBuckets, 0);

#pragma omp parallel default(none) shared(data, local_counts, n, chunk_size, exp, offset)
  {
    const int thread_id = omp_get_thread_num();
    const int start = thread_id * chunk_size;
    const int end = (start + chunk_size < n) ? (start + chunk_size) : n;
    int *local_count = &local_counts[thread_id * kBuckets];

    for (int i = start; i < end; i++) {
      int digit = ((data[i] + offset) / exp) % kBuckets;
      local_count[digit]++;
    }
  }

  std::array<int, kBuckets> global_count{};
  global_count.fill(0);

  for (int t = 0; t < num_threads; t++) {
    for (int b = 0; b < kBuckets; b++) {
      global_count[b] += local_counts[t * kBuckets + b];
    }
  }

  std::array<int, kBuckets> start_pos{};
  start_pos[0] = 0;
  for (int i = 1; i < kBuckets; i++) {
    start_pos[i] = start_pos[i - 1] + global_count[i - 1];
  }

  OutType output(n);

#pragma omp parallel default(none) shared(data, output, local_counts, n, chunk_size, start_pos, exp, offset)
  {
    const int thread_id = omp_get_thread_num();
    const int start = thread_id * chunk_size;
    const int end = (start + chunk_size < n) ? (start + chunk_size) : n;

    std::array<int, kBuckets> local_pos{};
    for (int b = 0; b < kBuckets; b++) {
      int base = start_pos[b];
      for (int t2 = 0; t2 < thread_id; t2++) {
        base += local_counts[t2 * kBuckets + b];
      }
      local_pos[b] = base;
    }

    for (int i = start; i < end; i++) {
      int digit = ((data[i] + offset) / exp) % kBuckets;
      output[local_pos[digit]++] = data[i];
    }
  }

  data = std::move(output);
}

void MelnikIRadixSortIntOMP::RadixSort(OutType &data) {
  if (data.empty()) {
    return;
  }

  int max_val = GetMaxValue(data);
  int min_val = *std::min_element(data.begin(), data.end());

  if (min_val >= 0) {
    for (int exp = 1; max_val / exp > 0; exp <<= kBitsPerDigit) {
      ParallelCountingSort(data, exp, 0);
    }
    return;
  }

  int offset = -min_val;
  for (int exp = 1; (max_val + offset) / exp > 0 || (min_val + offset) / exp > 0; exp <<= kBitsPerDigit) {
    ParallelCountingSort(data, exp, offset);
  }
}

}  // namespace melnik_i_radix_sort_int
