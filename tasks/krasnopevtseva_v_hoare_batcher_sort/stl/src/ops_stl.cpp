#include "krasnopevtseva_v_hoare_batcher_sort/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <future>
#include <stack>
#include <thread>
#include <utility>
#include <vector>

#include "krasnopevtseva_v_hoare_batcher_sort/common/include/common.hpp"

namespace krasnopevtseva_v_hoare_batcher_sort {

KrasnopevtsevaVHoareBatcherSortSTL::KrasnopevtsevaVHoareBatcherSortSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<int>();
}

bool KrasnopevtsevaVHoareBatcherSortSTL::ValidationImpl() {
  const auto &input = GetInput();
  return !input.empty();
}

bool KrasnopevtsevaVHoareBatcherSortSTL::PreProcessingImpl() {
  GetOutput() = std::vector<int>();
  return true;
}

bool KrasnopevtsevaVHoareBatcherSortSTL::RunImpl() {
  const auto &input = GetInput();
  std::size_t size = input.size();

  if (size <= 1) {
    GetOutput() = input;
    return true;
  }

  std::vector<int> res = input;

  int n = static_cast<int>(size);
  int numthreads = static_cast<int>(std::thread::hardware_concurrency());
  if (numthreads == 0) {
    numthreads = 1;
  }
  numthreads = std::min(n, numthreads);

  if (n < 1000) {
    QuickSort(res, 0, n - 1);
    GetOutput() = std::move(res);
    return true;
  }

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

  std::vector<std::thread> threads;
  threads.reserve(numthreads);

  for (int i = 0; i < numthreads; ++i) {
    threads.emplace_back([&res, &pointers, &sizes, i]() {
      int left = static_cast<int>(pointers[i] - res.data());
      int right = left + sizes[i] - 1;
      QuickSort(res, left, right);
    });
  }

  for (auto &t : threads) {
    t.join();
  }

  std::vector<int> temp(n);
  int step = 1;

  while (step < numthreads) {
    for (int i = 0; i < numthreads; i += 2 * step) {
      int left_block = i;
      int mid_block = std::min(i + step, numthreads);
      int right_block = std::min(i + 2 * step, numthreads);

      int left_start = (left_block == 0) ? 0 : static_cast<int>(pointers[left_block] - res.data());
      int left_end = static_cast<int>(pointers[mid_block] - res.data()) - 1;
      int right_start = left_end + 1;
      int right_end = (right_block == numthreads) ? n - 1 : static_cast<int>(pointers[right_block] - res.data()) - 1;

      int i1 = left_start, i2 = right_start, pos = left_start;
      while (i1 <= left_end && i2 <= right_end) {
        if (res[i1] <= res[i2]) {
          temp[pos++] = res[i1++];
        } else {
          temp[pos++] = res[i2++];
        }
      }
      while (i1 <= left_end) {
        temp[pos++] = res[i1++];
      }
      while (i2 <= right_end) {
        temp[pos++] = res[i2++];
      }
    }
    std::copy(temp.begin(), temp.end(), res.begin());
    step *= 2;
  }

  GetOutput() = std::move(res);
  return true;
}

bool KrasnopevtsevaVHoareBatcherSortSTL::PostProcessingImpl() {
  return true;
}

int KrasnopevtsevaVHoareBatcherSortSTL::Partition(std::vector<int> &arr, int first, int last) {
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

void KrasnopevtsevaVHoareBatcherSortSTL::InsertionSort(std::vector<int> &arr, int first, int last) {
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

void KrasnopevtsevaVHoareBatcherSortSTL::QuickSort(std::vector<int> &arr, int first, int last) {
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

void KrasnopevtsevaVHoareBatcherSortSTL::BatcherMergeBlocksStep(int *left_pointer, int &left_size, int *right_pointer,
                                                                int &right_size) {
  std::inplace_merge(left_pointer, right_pointer, right_pointer + right_size);
  left_size += right_size;
}

void KrasnopevtsevaVHoareBatcherSortSTL::BatcherMerge(int thread_input_size, std::vector<int *> &pointers,
                                                      std::vector<int> &sizes, int par_if_greater) {
  int pack = static_cast<int>(pointers.size());
  for (int step = 1; pack > 1; step *= 2, pack /= 2) {
    if ((thread_input_size / step) > par_if_greater) {
      std::vector<std::future<void>> futures;
      for (int off = 0; off < pack / 2; ++off) {
        futures.push_back(std::async(std::launch::async, [&pointers, &sizes, step, off]() {
          auto idx1 = static_cast<std::size_t>(2 * step) * static_cast<std::size_t>(off);
          auto idx2 = idx1 + static_cast<std::size_t>(step);
          BatcherMergeBlocksStep(pointers[idx1], sizes[idx1], pointers[idx2], sizes[idx2]);
        }));
      }
      for (auto &f : futures) {
        f.wait();
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

}  // namespace krasnopevtseva_v_hoare_batcher_sort
