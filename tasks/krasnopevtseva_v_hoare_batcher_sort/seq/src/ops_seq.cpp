#include "krasnopevtseva_v_hoare_batcher_sort/seq/include/ops_seq.hpp"

#include <cstddef>
#include <vector>

#include "krasnopevtseva_v_hoare_batcher_sort/common/include/common.hpp"
#include "util/include/util.hpp"

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
    quickBatcherSort(sort_v, 0, size - 1);
  }
  GetOutput() = sort_v;
  return true;
}

void KrasnopevtsevaVHoareBatcherSortSEQ::compareAndSwap(int &a, int &b) {
  if (a > b) {
    std::swap(a, b);
  }
}

void KrasnopevtsevaVHoareBatcherSortSEQ::batcherMerge(std::vector<int> &arr, int left, int right) {
  int n = right - left + 1;
  if (n <= 1) {
    return;
  }

  std::vector<int> temp(arr.begin() + left, arr.begin() + right + 1);

  std::function<void(int, int)> oddEvenMerge = [&](int l, int r) {
    if (l == r) {
      return;
    }

    int m = l + (r - l) / 2;
    oddEvenMerge(l, m);
    oddEvenMerge(m + 1, r);

    for (int i = l + 1; i + (m - l + 1) <= r; i += 2) {
      compareAndSwap(temp[i], temp[i + (m - l + 1)]);
    }
  };

  oddEvenMerge(0, n - 1);

  for (int i = 1; i + 1 < n; i += 2) {
    compareAndSwap(temp[i], temp[i + 1]);
  }
  for (int i = 0; i < n; i++) {
    arr[left + i] = temp[i];
  }
}

void KrasnopevtsevaVHoareBatcherSortSEQ::quickBatcherSort(std::vector<int> &arr, int left, int right) {
  if (left >= right) {
    return;
  }

  int mid = left + (right - left) / 2;
  if (arr[left] > arr[mid]) {
    std::swap(arr[left], arr[mid]);
  }
  if (arr[left] > arr[right]) {
    std::swap(arr[left], arr[right]);
  }
  if (arr[mid] > arr[right]) {
    std::swap(arr[mid], arr[right]);
  }

  std::swap(arr[mid], arr[right]);
  int pivot = arr[right];

  int i = left - 1;
  int j = right;

  while (true) {
    while (arr[++i] < pivot);
    while (arr[--j] > pivot);
    if (i >= j) {
      break;
    }
    std::swap(arr[i], arr[j]);
  }
  std::swap(arr[i], arr[right]);

  quickBatcherSort(arr, left, i - 1);
  quickBatcherSort(arr, i + 1, right);

  if (right - left > 32) {
    batcherMerge(arr, left, right);
  }
}
bool KrasnopevtsevaVHoareBatcherSortSEQ::PostProcessingImpl() {
  return true;
}
}  // namespace krasnopevtseva_v_hoare_batcher_sort
