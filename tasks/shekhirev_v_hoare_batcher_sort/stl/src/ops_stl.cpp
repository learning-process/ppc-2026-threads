#include "shekhirev_v_hoare_batcher_sort/stl/include/ops_stl.hpp"

#include <algorithm>
#include <limits>
#include <thread>
#include <utility>
#include <vector>

namespace shekhirev_v_hoare_batcher_sort {

namespace {

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
      int pivot = arr[l + (r - l) / 2];
      int i = l;
      int j = r;

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
  std::vector<int> buffer(chunk_size * 2);
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

  for (int x = 0; x < chunk_size; ++x) {
    arr[start1 + x] = buffer[x];
    arr[start2 + x] = buffer[chunk_size + x];
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
        [&data, i, chunk_size]() { OptimizedHoareSort(data, i * chunk_size, (i + 1) * chunk_size - 1); });
  }

  for (auto &t : threads) {
    t.join();
  }
  threads.clear();

  for (int p = 1; p < num_threads; p *= 2) {
    for (int k = p; k > 0; k /= 2) {
      for (int i = 0; i < num_threads - k; ++i) {
        if ((i / (p * 2)) == ((i + k) / (p * 2))) {
          int start_a = i * chunk_size;
          int start_b = (i + k) * chunk_size;

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

  if (padding > 0) {
    data.resize(orig_size);
  }

  return true;
}

bool ShekhirevHoareBatcherSortSTL::PostProcessingImpl() {
  return true;
}

}  // namespace shekhirev_v_hoare_batcher_sort
