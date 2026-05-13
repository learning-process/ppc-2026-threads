#include "shekhirev_v_hoare_batcher_sort/stl/include/ops_stl.hpp"

#include <algorithm>
#include <limits>
#include <thread>
#include <utility>
#include <vector>

#include "shekhirev_v_hoare_batcher_sort/common/include/common.hpp"

namespace shekhirev_v_hoare_batcher_sort {

namespace {

void SplitPartition(std::vector<int> &arr, int &l, int &r, int &i, int &j) {
  int pivot = arr[l + ((r - l) / 2)];
  i = l;
  j = r;

  while (i <= j) {
    while (arr[i] < pivot) {
      i++;
    }
    while (arr[j] > pivot) {
      j--;
    }
    if (i <= j) {
      std::swap(arr[i], arr[j]);
      i++;
      j--;
    }
  }
}

void OptimizedHoareSort(std::vector<int> &arr, int left, int right) {
  if (left >= right) {
    return;
  }

  std::vector<std::pair<int, int>> stack;
  stack.reserve(64);
  stack.emplace_back(left, right);

  while (!stack.empty()) {
    auto [l, r] = stack.back();
    stack.pop_back();

    while (l < r) {
      int i = 0;
      int j = 0;
      SplitPartition(arr, l, r, i, j);

      if (j - l < r - i) {
        if (i < r) {
          stack.emplace_back(i, r);
        }
        r = j;
      } else {
        if (l < j) {
          stack.emplace_back(l, j);
        }
        l = i;
      }
    }
  }
}

void MergeBlocks(std::vector<int> &arr, int start1, int start2, int chunk_size) {
  std::vector<int> buffer(static_cast<size_t>(chunk_size) * 2);
  int i = start1;
  int j = start2;
  int k = 0;

  int end1 = start1 + chunk_size;
  int end2 = start2 + chunk_size;

  while (i < end1 && j < end2) {
    if (arr[i] <= arr[j]) {
      buffer[k++] = arr[i++];
    } else {
      buffer[k++] = arr[j++];
    }
  }

  while (i < end1) {
    buffer[k++] = arr[i++];
  }
  while (j < end2) {
    buffer[k++] = arr[j++];
  }

  for (int idx = 0; idx < chunk_size; ++idx) {
    arr[start1 + idx] = buffer[idx];
    arr[start2 + idx] = buffer[chunk_size + idx];
  }
}

void BatcherMergePhase(std::vector<int> &data, int num_threads, int chunk_size) {
  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  for (int step_p = 1; step_p < num_threads; step_p *= 2) {
    for (int step_k = step_p; step_k > 0; step_k /= 2) {
      for (int idx = 0; idx < num_threads - step_k; ++idx) {
        if ((idx / (step_p * 2)) == ((idx + step_k) / (step_p * 2))) {
          int start_a = idx * chunk_size;
          int start_b = (idx + step_k) * chunk_size;

          threads.emplace_back(
              [&data, start_a, start_b, chunk_size]() { MergeBlocks(data, start_a, start_b, chunk_size); });
        }
      }

      for (auto &t : threads) {
        t.join();
      }
      threads.clear();
    }
  }
}

}  // namespace

ShekhirevHoareBatcherSortSTL::ShekhirevHoareBatcherSortSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ShekhirevHoareBatcherSortSTL::ValidationImpl() {
  return true;
}

bool ShekhirevHoareBatcherSortSTL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool ShekhirevHoareBatcherSortSTL::RunImpl() {
  auto &data = GetOutput();
  int orig_size = static_cast<int>(data.size());

  if (orig_size <= 1) {
    return true;
  }

  int hw_concurrency = static_cast<int>(std::thread::hardware_concurrency());
  if (hw_concurrency == 0) {
    hw_concurrency = 4;
  }

  int num_threads = 1;
  while (num_threads * 2 <= hw_concurrency && num_threads * 2 <= orig_size) {
    num_threads *= 2;
  }

  if (num_threads == 1) {
    OptimizedHoareSort(data, 0, orig_size - 1);
    return true;
  }

  int padding = (num_threads - (orig_size % num_threads)) % num_threads;
  if (padding > 0) {
    data.insert(data.end(), padding, std::numeric_limits<int>::max());
  }

  int total_size = static_cast<int>(data.size());
  int chunk_size = total_size / num_threads;

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  for (int i = 0; i < num_threads; ++i) {
    threads.emplace_back(
        [&data, i, chunk_size]() { OptimizedHoareSort(data, i * chunk_size, ((i + 1) * chunk_size) - 1); });
  }

  for (auto &t : threads) {
    t.join();
  }
  threads.clear();

  BatcherMergePhase(data, num_threads, chunk_size);

  if (padding > 0) {
    data.resize(orig_size);
  }

  return true;
}

bool ShekhirevHoareBatcherSortSTL::PostProcessingImpl() {
  return true;
}

}  // namespace shekhirev_v_hoare_batcher_sort
