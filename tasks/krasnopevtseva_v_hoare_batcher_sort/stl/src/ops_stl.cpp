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
  return !GetInput().empty();
}

bool KrasnopevtsevaVHoareBatcherSortSTL::PreProcessingImpl() {
  GetOutput() = std::vector<int>();
  return true;
}

int KrasnopevtsevaVHoareBatcherSortSTL::GetNumThreads(int n) {
  int numthreads = static_cast<int>(std::thread::hardware_concurrency());
  if (numthreads == 0) {
    numthreads = 1;
  }
  return std::min(n, numthreads);
}

void KrasnopevtsevaVHoareBatcherSortSTL::SetupChunks(std::vector<Chunk> &chunks, int n, int numthreads) {
  int thread_input_size = n / numthreads;
  int remainder = n % numthreads;

  int current_start = 0;
  for (int i = 0; i < numthreads; ++i) {
    int chunk_size = thread_input_size;
    if (i == numthreads - 1) {
      chunk_size += remainder;
    }
    chunks.push_back({nullptr, chunk_size, current_start, current_start + chunk_size - 1});
    current_start += chunk_size;
  }
}

void KrasnopevtsevaVHoareBatcherSortSTL::ParallelSortChunks(std::vector<int> &arr, std::vector<Chunk> &chunks) {
  std::vector<std::future<void>> futures;
  futures.reserve(chunks.size());

  for (auto &chunk : chunks) {
    futures.push_back(std::async(std::launch::async, [&arr, &chunk]() { QuickSort(arr, chunk.left, chunk.right); }));
  }

  for (auto &f : futures) {
    f.wait();
  }
}

void KrasnopevtsevaVHoareBatcherSortSTL::BatcherMergeBlocksStep(int *left_ptr, int &left_sz, int *right_ptr,
                                                                int &right_sz) {
  std::inplace_merge(left_ptr, right_ptr, right_ptr + right_sz);
  left_sz += right_sz;
}

void KrasnopevtsevaVHoareBatcherSortSTL::BatcherMergeLevel(int step, std::vector<Chunk> &chunks, int thread_input_size,
                                                           int par_if_greater) {
  int pack = static_cast<int>(chunks.size());
  bool do_parallel = (thread_input_size / step) > par_if_greater;

  auto merge_pair = [&](int off) {
    size_t idx1 = static_cast<size_t>(2 * step) * static_cast<size_t>(off);
    size_t idx2 = idx1 + static_cast<size_t>(step);
    BatcherMergeBlocksStep(chunks[idx1].ptr, chunks[idx1].size, chunks[idx2].ptr, chunks[idx2].size);
    chunks[idx1].right = chunks[idx2].right;
  };

  if (do_parallel) {
    std::vector<std::future<void>> futures;
    futures.reserve(pack / 2);
    for (int off = 0; off < pack / 2; ++off) {
      futures.push_back(std::async(std::launch::async, merge_pair, off));
    }
    for (auto &f : futures) {
      f.wait();
    }
  } else {
    for (int off = 0; off < pack / 2; ++off) {
      merge_pair(off);
    }
  }
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
  int numthreads = GetNumThreads(n);

  std::vector<Chunk> chunks;
  SetupChunks(chunks, n, numthreads);

  for (auto &chunk : chunks) {
    chunk.ptr = res.data() + chunk.left;
  }

  ParallelSortChunks(res, chunks);

  for (int step = 1; static_cast<int>(chunks.size()) > 1; step *= 2) {
    int pack = static_cast<int>(chunks.size());
    BatcherMergeLevel(step, chunks, n / numthreads, 32);

    if ((pack / 2) - 1 == 0) {
      BatcherMergeBlocksStep(chunks[0].ptr, chunks[0].size, chunks.back().ptr, chunks.back().size);
      chunks[0].right = chunks.back().right;
      chunks.resize(1);
    } else if ((pack / 2) % 2 != 0) {
      size_t idx1 = static_cast<size_t>(2 * step) * static_cast<size_t>((pack / 2) - 2);
      size_t idx2 = idx1 + static_cast<size_t>(step);
      BatcherMergeBlocksStep(chunks[idx1].ptr, chunks[idx1].size, chunks[idx2].ptr, chunks[idx2].size);
      chunks.erase(chunks.begin() + static_cast<ptrdiff_t>(idx2));
    } else {
      int new_size = pack / 2;
      chunks.resize(new_size);
    }
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

}  // namespace krasnopevtseva_v_hoare_batcher_sort
