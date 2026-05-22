#include "krasnopevtseva_v_hoare_batcher_sort/all/include/ops_all.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stack>
#include <utility>
#include <vector>

#include "krasnopevtseva_v_hoare_batcher_sort/common/include/common.hpp"

namespace krasnopevtseva_v_hoare_batcher_sort {

KrasnopevtsevaVHoareBatcherSortALL::KrasnopevtsevaVHoareBatcherSortALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<int>();
}

bool KrasnopevtsevaVHoareBatcherSortALL::ValidationImpl() {
  const auto &input = GetInput();
  return !input.empty();
}

bool KrasnopevtsevaVHoareBatcherSortALL::PreProcessingImpl() {
  GetOutput() = std::vector<int>();
  return true;
}

bool KrasnopevtsevaVHoareBatcherSortALL::PostProcessingImpl() {
  return true;
}

int KrasnopevtsevaVHoareBatcherSortALL::Partition(std::vector<int> &arr, int first, int last) {
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

void KrasnopevtsevaVHoareBatcherSortALL::InsertionSort(std::vector<int> &arr, int first, int last) {
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

void KrasnopevtsevaVHoareBatcherSortALL::QuickSort(std::vector<int> &arr, int first, int last) {
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

void KrasnopevtsevaVHoareBatcherSortALL::BatcherMergeBlocksStep(int *left_pointer, int &left_size, int *right_pointer,
                                                                int &right_size) {
  std::inplace_merge(left_pointer, right_pointer, right_pointer + right_size);
  left_size += right_size;
}

void KrasnopevtsevaVHoareBatcherSortALL::BatcherMerge(int thread_input_size, std::vector<int *> &pointers,
                                                      std::vector<int> &sizes, int par_if_greater) {
  int pack = static_cast<int>(pointers.size());
  for (int step = 1; pack > 1; step *= 2, pack /= 2) {
    bool do_parallel = (thread_input_size / step) > par_if_greater;

    if (do_parallel) {
#pragma omp parallel for default(none) shared(pointers, sizes, pack, step)
      for (int off = 0; off < pack / 2; ++off) {
        auto idx1 = static_cast<std::size_t>(2 * step) * static_cast<std::size_t>(off);
        auto idx2 = idx1 + static_cast<std::size_t>(step);
        BatcherMergeBlocksStep(pointers[idx1], sizes[idx1], pointers[idx2], sizes[idx2]);
      }
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

void KrasnopevtsevaVHoareBatcherSortALL::ParallelSortChunks(std::vector<int> &arr, int n, int numthreads) {
  if (n <= 0) {
    return;
  }

  numthreads = std::min(n, numthreads);
  if (numthreads <= 0) {
    numthreads = 1;
  }

  int thread_input_size = n / numthreads;
  int thread_input_remainder_size = n % numthreads;

  std::vector<int *> pointers(numthreads);
  std::vector<int> sizes(numthreads);

  for (int i = 0; i < numthreads; ++i) {
    std::ptrdiff_t offset = static_cast<std::ptrdiff_t>(i) * static_cast<std::ptrdiff_t>(thread_input_size);
    pointers[i] = arr.data() + offset;
    sizes[i] = thread_input_size;
  }
  sizes.back() += thread_input_remainder_size;

#pragma omp parallel for default(none) shared(arr, pointers, sizes, numthreads)
  for (int i = 0; i < numthreads; ++i) {
    int left = static_cast<int>(pointers[i] - arr.data());
    int right = left + sizes[i] - 1;
    if (left < right) {
      QuickSort(arr, left, right);
    }
  }

  BatcherMerge(thread_input_size, pointers, sizes, 32);
}

void KrasnopevtsevaVHoareBatcherSortALL::SortLocalData(std::vector<int> &data) {
  int n = static_cast<int>(data.size());
  if (n <= 0) {
    return;
  }

  int numthreads = omp_get_max_threads();
  numthreads = std::min(n, numthreads);

  if (n < 1000) {
    QuickSort(data, 0, n - 1);
  } else {
    ParallelSortChunks(data, n, numthreads);
  }
}

bool KrasnopevtsevaVHoareBatcherSortALL::RunImpl() {
  const auto &input = GetInput();
  std::vector<int> result = input;

  if (result.size() <= 1) {
    GetOutput() = result;
    return true;
  }

  SortLocalData(result);
  GetOutput() = std::move(result);
  return true;
}

}  // namespace krasnopevtseva_v_hoare_batcher_sort
