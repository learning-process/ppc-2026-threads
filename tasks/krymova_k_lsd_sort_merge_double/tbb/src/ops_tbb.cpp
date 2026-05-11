#include "krymova_k_lsd_sort_merge_double/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <vector>

#include "krymova_k_lsd_sort_merge_double/common/include/common.hpp"

namespace krymova_k_lsd_sort_merge_double {

KrymovaKLsdSortMergeDoubleTBB::KrymovaKLsdSortMergeDoubleTBB(const InType &in)
    : num_threads_(tbb::this_task_arena::max_concurrency()) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool KrymovaKLsdSortMergeDoubleTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool KrymovaKLsdSortMergeDoubleTBB::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

uint64_t KrymovaKLsdSortMergeDoubleTBB::DoubleToULL(double d) {
  uint64_t ull = 0;
  std::memcpy(&ull, &d, sizeof(double));

  if ((ull & 0x8000000000000000ULL) != 0U) {
    ull = ~ull;
  } else {
    ull |= 0x8000000000000000ULL;
  }

  return ull;
}

double KrymovaKLsdSortMergeDoubleTBB::ULLToDouble(uint64_t ull) {
  if ((ull & 0x8000000000000000ULL) != 0U) {
    ull &= 0x7FFFFFFFFFFFFFFFULL;
  } else {
    ull = ~ull;
  }

  double d = 0.0;
  std::memcpy(&d, &ull, sizeof(double));
  return d;
}

void KrymovaKLsdSortMergeDoubleTBB::LSDSortDouble(double *arr, int size) {
  if (size <= 1) {
    return;
  }

  const int k_bits_per_pass = 8;
  const int k_radix = 1 << k_bits_per_pass;
  const int k_passes = static_cast<int>(sizeof(double)) * 8 / k_bits_per_pass;

  std::vector<uint64_t> ull_arr(size);
  std::vector<uint64_t> ull_tmp(size);

  for (int i = 0; i < size; ++i) {
    ull_arr[i] = DoubleToULL(arr[i]);
  }

  std::vector<unsigned int> count(k_radix, 0U);

  for (int pass = 0; pass < k_passes; ++pass) {
    int shift = pass * k_bits_per_pass;

    std::ranges::fill(count, 0U);

    for (int i = 0; i < size; ++i) {
      unsigned int digit = (ull_arr[i] >> shift) & (k_radix - 1);
      ++count[digit];
    }

    for (int i = 1; i < k_radix; ++i) {
      count[i] += count[i - 1];
    }

    for (int i = size - 1; i >= 0; --i) {
      unsigned int digit = (ull_arr[i] >> shift) & (k_radix - 1);
      ull_tmp[--count[digit]] = ull_arr[i];
    }

    ull_arr.swap(ull_tmp);
  }

  for (int i = 0; i < size; ++i) {
    arr[i] = ULLToDouble(ull_arr[i]);
  }
}

void KrymovaKLsdSortMergeDoubleTBB::MergeSections(double *left, const double *right, int left_size, int right_size) {
  std::vector<double> temp(left_size);
  std::ranges::copy(left, left + left_size, temp.begin());

  int left_index = 0;
  int right_index = 0;
  int dest_index = 0;

  while (left_index < left_size && right_index < right_size) {
    if (temp[left_index] <= right[right_index]) {
      left[dest_index++] = temp[left_index++];
    } else {
      left[dest_index++] = right[right_index++];
    }
  }

  while (left_index < left_size) {
    left[dest_index++] = temp[left_index++];
  }
}

void KrymovaKLsdSortMergeDoubleTBB::SortSectionsParallel(double *arr, int size, int portion) {
  int num_blocks = (size + portion - 1) / portion;
  const size_t grain_size = 1;
  tbb::parallel_for(tbb::blocked_range<int>(0, num_blocks, grain_size), [&](const tbb::blocked_range<int> &range) {
    for (int block_index = range.begin(); block_index != range.end(); ++block_index) {
      int start = block_index * portion;
      int current_size = std::min(portion, size - start);
      LSDSortDouble(arr + start, current_size);
    }
  });
}

void KrymovaKLsdSortMergeDoubleTBB::IterativeMergeSort(double *arr, int size, int portion) {
  if (size <= 1) {
    return;
  }

  SortSectionsParallel(arr, size, portion);

  for (int merge_size = portion; merge_size < size; merge_size *= 2) {
    int num_pairs = (size + 2 * merge_size - 1) / (2 * merge_size);
    const size_t grain_size = 16;
    tbb::parallel_for(tbb::blocked_range<int>(0, num_pairs, grain_size), [&](const tbb::blocked_range<int> &range) {
      for (int pair_index = range.begin(); pair_index != range.end(); ++pair_index) {
        int start = pair_index * 2 * merge_size;
        int left_size = merge_size;
        int right_size = std::min(merge_size, size - (start + merge_size));
        if (right_size > 0) {
          MergeSections(arr + start, arr + start + left_size, left_size, right_size);
        }
      }
    });
  }
}

bool KrymovaKLsdSortMergeDoubleTBB::RunImpl() {
  OutType &output = GetOutput();
  int size = static_cast<int>(output.size());

  if (size <= 1) {
    return true;
  }

  int portion = std::max(1, size / num_threads_);
  portion = std::min(portion, 1000000);

  IterativeMergeSort(output.data(), size, portion);

  return true;
}

bool KrymovaKLsdSortMergeDoubleTBB::PostProcessingImpl() {
  const OutType &output = GetOutput();

  for (size_t i = 1; i < output.size(); ++i) {
    if (output[i] < output[i - 1]) {
      return false;
    }
  }

  return true;
}

}  // namespace krymova_k_lsd_sort_merge_double
