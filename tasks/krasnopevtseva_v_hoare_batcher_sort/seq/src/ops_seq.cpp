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
  stack.push({left, right});

  while (!stack.empty()) {
    auto [l, r] = stack.top();
    stack.pop();

    if (l >= r) {
      continue;
    }
    if (r - l < 16) {
      for (int i = l + 1; i <= r; i++) {
        int key = arr[i];
        int j = i - 1;
        while (j >= l && arr[j] > key) {
          arr[j + 1] = arr[j];
          j--;
        }
        arr[j + 1] = key;
      }
      continue;
    }

    int mid = l + (r - l) / 2;
    if (arr[l] > arr[mid]) {
      std::swap(arr[l], arr[mid]);
    }
    if (arr[l] > arr[r]) {
      std::swap(arr[l], arr[r]);
    }
    if (arr[mid] > arr[r]) {
      std::swap(arr[mid], arr[r]);
    }

    std::swap(arr[mid], arr[r - 1]);
    int pivot = arr[r - 1];

    int i = l;
    int j = r - 1;

    while (true) {
      while (arr[++i] < pivot);
      while (arr[--j] > pivot);
      if (i >= j) {
        break;
      }
      std::swap(arr[i], arr[j]);
    }

    std::swap(arr[i], arr[r - 1]);

    if (i - l < r - i) {
      stack.push({i + 1, r});
      stack.push({l, i - 1});
    } else {
      stack.push({l, i - 1});
      stack.push({i + 1, r});
    }
  }

  if (right - left > 32) {
    BatcherMerge(arr, left, right);
  }
}
bool KrasnopevtsevaVHoareBatcherSortSEQ::PostProcessingImpl() {
  return true;
}
}  // namespace krasnopevtseva_v_hoare_batcher_sort
