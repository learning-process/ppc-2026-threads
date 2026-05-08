#include "krasnopevtseva_v_hoare_batcher_sort/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cstddef>
#include <stack>
#include <thread>
#include <utility>
#include <vector>

#include "krasnopevtseva_v_hoare_batcher_sort/common/include/common.hpp"
#include "oneapi/tbb/parallel_for.h"

namespace krasnopevtseva_v_hoare_batcher_sort {

KrasnopevtsevaVHoareBatcherSortTBB::KrasnopevtsevaVHoareBatcherSortTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<int>();
}

bool KrasnopevtsevaVHoareBatcherSortTBB::ValidationImpl() {
  const auto &input = GetInput();
  return !input.empty();
}

bool KrasnopevtsevaVHoareBatcherSortTBB::PreProcessingImpl() {
  GetOutput() = std::vector<int>();
  return true;
}

bool KrasnopevtsevaVHoareBatcherSortTBB::RunImpl() {
  const auto &input = GetInput();
  std::size_t size = input.size();

  if (size <= 1) {
    GetOutput() = input;
    return true;
  }

  int numthreads = static_cast<int>(std::thread::hardware_concurrency());
  static tbb::global_control control(tbb::global_control::max_allowed_parallelism, numthreads);

  std::vector<int> res = input;

  int n = static_cast<int>(size);
  numthreads = std::min(n, numthreads);

  int thread_input_size = n / numthreads;
  int thread_input_remainder_size = n % numthreads;

  std::vector<int *> pointers(numthreads);
  std::vector<int> sizes(numthreads);
  for (int i = 0; i < numthreads; ++i) {
    auto offset = static_cast<std::ptrdiff_t>(i) * static_cast<std::ptrdiff_t>(thread_input_size);
    pointers[i] = res.data() + offset;
    sizes[i] = thread_input_size;
  }
  sizes[sizes.size() - 1] += thread_input_remainder_size;

  tbb::parallel_for(tbb::blocked_range<int>(0, numthreads, 1), [&](const tbb::blocked_range<int> &r) {
    for (int i = r.begin(); i < r.end(); ++i) {
      int left = static_cast<int>(pointers[i] - res.data());
      int right = left + sizes[i] - 1;
      QuickSort(res, left, right);
    }
  }, tbb::simple_partitioner());

  BatcherMerge(thread_input_size, pointers, sizes, 32);

  GetOutput() = std::move(res);
  return true;
}

bool KrasnopevtsevaVHoareBatcherSortTBB::PostProcessingImpl() {
  return true;
}

int KrasnopevtsevaVHoareBatcherSortTBB::Partition(std::vector<int> &arr, int first, int last) {
  int i = first - 1;
  int value = arr[last];

  for (int j = first; j <= last - 1; ++j) {
    if (arr[j] <= value) {
      ++i;
      std::swap(arr[i], arr[j]);
    }
  }
  std::swap(arr[i + 1], arr[last]);
  return i + 1;
}

void KrasnopevtsevaVHoareBatcherSortTBB::InsertionSort(std::vector<int> &arr, int first, int last) {
  for (int i = first + 1; i <= last; ++i) {
    int key = arr[i];
    int j = i - 1;
    while (j >= first && arr[j] > key) {
      arr[j + 1] = arr[j];
      --j;
    }
    arr[j + 1] = key;
  }
}

void KrasnopevtsevaVHoareBatcherSortTBB::QuickSort(std::vector<int> &arr, int first, int last) {
  std::stack<std::pair<int, int>> stack;
  stack.emplace(first, last);

  while (!stack.empty()) {
    auto [l, r] = stack.top();
    stack.pop();

    if (l >= r) {
      continue;
    }

    if (r - l < 16) {
      InsertionSort(arr, l, r);
      continue;
    }

    int iter = Partition(arr, l, r);

    if (iter - l < r - iter) {
      stack.emplace(iter + 1, r);
      stack.emplace(l, iter - 1);
    } else {
      stack.emplace(l, iter - 1);
      stack.emplace(iter + 1, r);
    }
  }
}

void KrasnopevtsevaVHoareBatcherSortTBB::BatcherMergeBlocksStep(int *left_pointer, int &left_size, int *right_pointer,
                                                                int &right_size) {
  std::inplace_merge(left_pointer, right_pointer, right_pointer + right_size);
  left_size += right_size;
}

void KrasnopevtsevaVHoareBatcherSortTBB::BatcherMerge(int thread_input_size, std::vector<int *> &pointers,
                                                      std::vector<int> &sizes, int par_if_greater) {
  int pack = static_cast<int>(pointers.size());
  for (int step = 1; pack > 1; step *= 2, pack /= 2) {
    if ((thread_input_size / step) > par_if_greater) {
      tbb::parallel_for(tbb::blocked_range<int>(0, pack / 2, 1), [&](const tbb::blocked_range<int> &r) {
        for (int off = r.begin(); off < r.end(); ++off) {
          auto idx1 = static_cast<std::size_t>(2 * step) * static_cast<std::size_t>(off);
          auto idx2 = idx1 + static_cast<std::size_t>(step);
          BatcherMergeBlocksStep(pointers[idx1], sizes[idx1], pointers[idx2], sizes[idx2]);
        }
      }, tbb::simple_partitioner());
    } else {
      for (int off = 0; off < pack / 2; ++off) {
        auto idx1 = static_cast<std::size_t>(2 * step) * static_cast<std::size_t>(off);
        auto idx2 = idx1 + static_cast<std::size_t>(step);
        BatcherMergeBlocksStep(pointers[idx1], sizes[idx1], pointers[idx2], sizes[idx2]);
      }
    }

    if ((pack / 2) - 1 == 0) {
      BatcherMergeBlocksStep(pointers[0], sizes[sizes.size() - 1], pointers[pointers.size() - 1],
                             sizes[sizes.size() - 1]);
    } else if ((pack / 2) % 2 != 0) {
      auto idx1 = static_cast<std::size_t>(2 * step) * static_cast<std::size_t>((pack / 2) - 2);
      auto idx2 = static_cast<std::size_t>(2 * step) * static_cast<std::size_t>((pack / 2) - 1);
      BatcherMergeBlocksStep(pointers[idx1], sizes[idx1], pointers[idx2], sizes[idx2]);
    }
  }
}

}  // namespace krasnopevtseva_v_hoare_batcher_sort
