#include "krasnopevtseva_v_hoare_batcher_sort/seq/include/ops_seq.hpp"

#include <cstddef>
#include <functional>
#include <stack>
#include <utility>
#include <vector>

#include "krasnopevtseva_v_hoare_batcher_sort/common/include/common.hpp"

namespace krasnopevtseva_v_hoare_batcher_sort {

KrasnopevtsevaVHoareBatcherSortSEQ::KrasnopevtsevaVHoareBatcherSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<int>();
}

bool KrasnopevtsevaVHoareBatcherSortSEQ::ValidationImpl() {
  const auto &input = GetInput();
  return (!input.empty());
}

bool KrasnopevtsevaVHoareBatcherSortSEQ::PreProcessingImpl() {
  GetOutput() = std::vector<int>();
  return true;
}

bool KrasnopevtsevaVHoareBatcherSortSEQ::RunImpl() {
  const auto &input = GetInput();
  size_t size = input.size();
  std::vector<int> sort_v = input;

  if (size > 1) {
    QuickBatcherSort(sort_v, 0, static_cast<int>(size - 1));
  }
  GetOutput() = sort_v;
  return true;
}

void KrasnopevtsevaVHoareBatcherSortSEQ::CompareAndSwap(int &a, int &b) {
  if (a > b) {
    std::swap(a, b);
  }
}

void KrasnopevtsevaVHoareBatcherSortSEQ::BatcherMerge(std::vector<int> &arr, int left, int right) {
  int n = right - left + 1;
  if (n <= 1) {
    return;
  }

  std::vector<int> temp(arr.begin() + left, arr.begin() + right + 1);

  std::function<void(int, int)> odd_even_merge = [&](int l, int r) {
    if (l == r) {
      return;
    }

    int m = l + ((r - l) / 2);
    odd_even_merge(l, m);
    odd_even_merge(m + 1, r);

    for (int i = l + 1; i + (m - l + 1) <= r; i += 2) {
      CompareAndSwap(temp[i], temp[i + (m - l + 1)]);
    }
  };

  odd_even_merge(0, n - 1);

  for (int i = 1; i + 1 < n; i += 2) {
    CompareAndSwap(temp[i], temp[i + 1]);
  }
  for (int i = 0; i < n; i++) {
    arr[left + i] = temp[i];
  }
}

void KrasnopevtsevaVHoareBatcherSortSEQ::QuickBatcherSort(std::vector<int> &arr, int left, int right) {
  std::stack<std::pair<int, int>> stack;
  stack.emplace(left, right);

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

    int pivot_index = Partition(arr, l, r);

    if (pivot_index - l < r - pivot_index) {
      stack.emplace(pivot_index + 1, r);
      stack.emplace(l, pivot_index - 1);
    } else {
      stack.emplace(l, pivot_index - 1);
      stack.emplace(pivot_index + 1, r);
    }
  }

  if (right - left > 32) {
    BatcherMerge(arr, left, right);
  }
}

void KrasnopevtsevaVHoareBatcherSortSEQ::InsertionSort(std::vector<int> &arr, int left, int right) {
  for (int i = left + 1; i <= right; ++i) {
    int key = arr[i];
    int j = i - 1;
    while (j >= left && arr[j] > key) {
      arr[j + 1] = arr[j];
      --j;
    }
    arr[j + 1] = key;
  }
}

int KrasnopevtsevaVHoareBatcherSortSEQ::Partition(std::vector<int> &arr, int left, int right) {
  int mid = left + ((right - left) / 2);
  if (arr[left] > arr[mid]) {
    std::swap(arr[left], arr[mid]);
  }
  if (arr[left] > arr[right]) {
    std::swap(arr[left], arr[right]);
  }
  if (arr[mid] > arr[right]) {
    std::swap(arr[mid], arr[right]);
  }

  std::swap(arr[mid], arr[right - 1]);
  int pivot = arr[right - 1];

  int i = left;
  int j = right - 1;

  while (true) {
    while (arr[++i] < pivot) {
    }
    while (arr[--j] > pivot) {
    }
    if (i >= j) {
      break;
    }
    std::swap(arr[i], arr[j]);
  }

  std::swap(arr[i], arr[right - 1]);
  return i;
}

bool KrasnopevtsevaVHoareBatcherSortSEQ::PostProcessingImpl() {
  return true;
}
}  // namespace krasnopevtseva_v_hoare_batcher_sort
