#include "krasnopevtseva_v_hoare_batcher_sort/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stack>
#include <utility>
#include <vector>

#include "krasnopevtseva_v_hoare_batcher_sort/common/include/common.hpp"
namespace krasnopevtseva_v_hoare_batcher_sort {
KrasnopevtsevaVHoareBatcherSortOMP::KrasnopevtsevaVHoareBatcherSortOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<int>();
}

bool KrasnopevtsevaVHoareBatcherSortOMP::ValidationImpl() {
  const auto &input = GetInput();
  return !input.empty();
}

bool KrasnopevtsevaVHoareBatcherSortOMP::PreProcessingImpl() {
  GetOutput() = std::vector<int>();
  return true;
}

bool KrasnopevtsevaVHoareBatcherSortOMP::RunImpl() {
  const auto &input = GetInput();
  std::size_t size = input.size();

  if (size <= 1) {
    GetOutput() = input;
    return true;
  }

  std::vector<int> res = input;

  int n = static_cast<int>(size);
  int numthreads = omp_get_max_threads();
  numthreads = std::min(n, numthreads);

  int thread_input_size = n / numthreads;
  int thread_input_remainder_size = n % numthreads;

  std::vector<int *> pointers(numthreads);
  std::vector<int> sizes(numthreads);
  for (int i = 0; i < numthreads; ++i) {
    std::ptrdiff_t offset = static_cast<std::ptrdiff_t>(i) * static_cast<std::ptrdiff_t>(thread_input_size);
    pointers[i] = res.data() + offset;
    sizes[i] = thread_input_size;
  }
  sizes[sizes.size() - 1] += thread_input_remainder_size;

#pragma omp parallel for default(none) shared(res, pointers, sizes, numthreads)
  for (int i = 0; i < numthreads; ++i) {
    int left = static_cast<int>(pointers[i] - res.data());
    int right = left + sizes[i] - 1;
    QuickSort(res, left, right);
  }

  BatcherMerge(thread_input_size, pointers, sizes, 32);

  GetOutput() = std::move(res);
  return true;
}

bool KrasnopevtsevaVHoareBatcherSortOMP::PostProcessingImpl() {
  return true;
}

int KrasnopevtsevaVHoareBatcherSortOMP::Partition(std::vector<int> &arr, int first, int last) {
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

void KrasnopevtsevaVHoareBatcherSortOMP::InsertionSort(std::vector<int> &arr, int first, int last) {
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

void KrasnopevtsevaVHoareBatcherSortOMP::QuickSort(std::vector<int> &arr, int first, int last) {
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

void KrasnopevtsevaVHoareBatcherSortOMP::BatcherMergeBlocksStep(int *left_pointer, int &left_size, int *right_pointer,
                                                                int &right_size) {
  std::inplace_merge(left_pointer, right_pointer, right_pointer + right_size);
  left_size += right_size;
}

void KrasnopevtsevaVHoareBatcherSortOMP::BatcherMerge(int thread_input_size, std::vector<int *> &pointers,
                                                      std::vector<int> &sizes, int par_if_greater) {
  int pack = static_cast<int>(pointers.size());
  for (int step = 1; pack > 1; step *= 2, pack /= 2) {
#pragma omp parallel for default(none) shared(pointers, sizes, pack, step, thread_input_size, \
                                                  par_if_greater) if ((thread_input_size / step) > par_if_greater)
    for (int off = 0; off < pack / 2; ++off) {
      auto idx1 = static_cast<std::size_t>(2 * step) * static_cast<std::size_t>(off);
      auto idx2 = (static_cast<std::size_t>(2 * step) * static_cast<std::size_t>(off)) + static_cast<std::size_t>(step);
      BatcherMergeBlocksStep(pointers[idx1], sizes[idx1], pointers[idx2], sizes[idx2]);
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
