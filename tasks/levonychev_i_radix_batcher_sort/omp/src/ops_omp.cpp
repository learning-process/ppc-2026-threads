#include "levonychev_i_radix_batcher_sort/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <vector>

#include "levonychev_i_radix_batcher_sort/common/include/common.hpp"
#include "util/include/util.hpp"

namespace levonychev_i_radix_batcher_sort {

LevonychevIRadixBatcherSortOMP::LevonychevIRadixBatcherSortOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

void LevonychevIRadixBatcherSortOMP::CountingSort(InType &arr, size_t byte_index) {
  const size_t byte = 256;
  std::vector<int> count(byte, 0);
  OutType result(arr.size());

  bool is_last_byte = (byte_index == (sizeof(int) - 1ULL));
  for (auto number : arr) {
    int value_of_byte = (number >> (byte_index * 8ULL)) & 0xFF;

    if (is_last_byte) {
      value_of_byte ^= 0x80;
    }

    ++count[value_of_byte];
  }

  for (size_t i = 1ULL; i < byte; ++i) {
    count[i] += count[i - 1];
  }

  for (int &val : std::ranges::reverse_view(arr)) {
    int value_of_byte = (val >> (byte_index * 8ULL)) & 0xFF;

    if (is_last_byte) {
      value_of_byte ^= 0x80;
    }

    result[--count[value_of_byte]] = val;
  }
  arr = result;
}

bool LevonychevIRadixBatcherSortOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool LevonychevIRadixBatcherSortOMP::PreProcessingImpl() {
  return true;
}

void LevonychevIRadixBatcherSortOMP::BatcherCompareRange(std::vector<int> &arr, int j, int k, int p2) {
  int range = std::min(k, static_cast<int>(arr.size()) - j - k);
  for (int i = 0; i < range; ++i) {
    int idx1 = j + i;
    int idx2 = j + i + k;
    if ((idx1 & p2) == (idx2 & p2) && (arr[idx1] > arr[idx2])) {
      std::swap(arr[idx1], arr[idx2]);
    }
  }
}

void LevonychevIRadixBatcherSortOMP::BatcherMergeIterative(std::vector<int> &arr, int start_p, int threads) {
  int n = static_cast<int>(arr.size());
  for (int pv = start_p; pv < n; pv <<= 1) {
    int p2 = pv << 1;
    for (int k = pv; k > 0; k >>= 1) {
#pragma omp parallel for schedule(static) default(none) shared(n, pv, p2, k, arr) num_threads(threads)
      for (int j = k % pv; j < n - k; j += 2 * k) {
        BatcherCompareRange(arr, j, k, p2);
      }
    }
  }
}

bool LevonychevIRadixBatcherSortOMP::RunImpl() {
  GetOutput() = GetInput();
  int n = static_cast<int>(GetOutput().size());
  if (n <= 1) {
    return true;
  }

  int num_threads = ppc::util::GetNumThreads();
  int block_size = n / num_threads;

#pragma omp parallel default(none) shared(num_threads, block_size, n, std::ranges::copy) num_threads(num_threads)
  {
    int tid = omp_get_thread_num();
    int left = tid * block_size;
    int right = (tid == num_threads - 1) ? n : (tid + 1) * block_size;

    if (left < right) {
      std::vector<int> local_block(GetOutput().begin() + left, GetOutput().begin() + right);
      for (size_t i = 0; i < sizeof(int); ++i) {
        CountingSort(local_block, i);
      }
      std::ranges::copy(local_block, GetOutput().begin() + left);
    }
  }

  int start_p = 1;
  while (start_p < block_size) {
    start_p <<= 1;
  }
  BatcherMergeIterative(GetOutput(), start_p, num_threads);

  return true;
}

bool LevonychevIRadixBatcherSortOMP::PostProcessingImpl() {
  return true;
}

}  // namespace levonychev_i_radix_batcher_sort
