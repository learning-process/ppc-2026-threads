#include "levonychev_i_radix_batcher_sort/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <vector>

#include "levonychev_i_radix_batcher_sort/common/include/common.hpp"
#include "util/include/util.hpp"

namespace levonychev_i_radix_batcher_sort {

LevonychevIRadixBatcherSortTBB::LevonychevIRadixBatcherSortTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

void LevonychevIRadixBatcherSortTBB::CountingSort(InType &arr, size_t byte_index) {
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

bool LevonychevIRadixBatcherSortTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool LevonychevIRadixBatcherSortTBB::PreProcessingImpl() {
  return true;
}

void LevonychevIRadixBatcherSortTBB::BatcherCompareRange(std::vector<int> &arr, int j, int k, int p2) {
  int range = std::min(k, static_cast<int>(arr.size()) - j - k);
  for (int i = 0; i < range; ++i) {
    int idx1 = j + i;
    int idx2 = j + i + k;
    if ((idx1 & p2) == (idx2 & p2) && (arr[idx1] > arr[idx2])) {
      std::swap(arr[idx1], arr[idx2]);
    }
  }
}

void LevonychevIRadixBatcherSortTBB::BatcherMergeIterative(std::vector<int> &arr, int start_p) {
  int n = static_cast<int>(arr.size());
  for (int pv = start_p; pv < n; pv <<= 1) {
    int p2 = pv << 1;
    for (int k = pv; k > 0; k >>= 1) {
      int num_iters = (n - k - (k % pv) + 2 * k - 1) / (2 * k);

      tbb::parallel_for(0, num_iters, [&](int i) {
        int j = (k % pv) + (i * (2 * k));
        BatcherCompareRange(arr, j, k, p2);
      });
    }
  }
}

bool LevonychevIRadixBatcherSortTBB::RunImpl() {
  GetOutput() = GetInput();
  int n = static_cast<int>(GetOutput().size());
  if (n <= 1) {
    return true;
  }
  int num_threads = ppc::util::GetNumThreads();
  int grain = std::max(1, n / num_threads);

  tbb::parallel_for(tbb::blocked_range<int>(0, n, grain), [&](const tbb::blocked_range<int> &r) {
    int left = r.begin();
    int right = r.end();

    std::vector<int> local_block(GetOutput().begin() + left, GetOutput().begin() + right);
    for (size_t i = 0; i < sizeof(int); ++i) {
      CountingSort(local_block, i);
    }
    std::ranges::copy(local_block, GetOutput().begin() + left);
  });

  int block_size = grain;
  int start_p = 1;
  while (start_p < block_size) {
    start_p <<= 1;
  }

  BatcherMergeIterative(GetOutput(), start_p);
  return true;
}

bool LevonychevIRadixBatcherSortTBB::PostProcessingImpl() {
  return true;
}

}  // namespace levonychev_i_radix_batcher_sort
