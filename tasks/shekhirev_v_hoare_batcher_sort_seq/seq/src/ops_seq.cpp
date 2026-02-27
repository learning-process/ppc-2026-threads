#include "seq/include/ops_seq.hpp"

#include <algorithm>
#include <climits>
#include <utility>

namespace shekhirev_v_hoare_batcher_sort_seq {

ShekhirevHoareBatcherSortSEQ::ShekhirevHoareBatcherSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ShekhirevHoareBatcherSortSEQ::ValidationImpl() {
  return true;
}

bool ShekhirevHoareBatcherSortSEQ::PreProcessingImpl() {
  const auto &in = GetInput();
  size_t original_size = in.size();

  size_t p2 = 1;
  while (p2 < original_size) {
    p2 *= 2;
  }

  GetOutput() = in;
  GetOutput().resize(p2, INT_MAX);
  return true;
}

void ShekhirevHoareBatcherSortSEQ::HoareSort(std::vector<int> &arr, int left, int right) {
  if (left >= right) {
    return;
  }
  int pivot = arr[left + (right - left) / 2];
  int i = left;
  int j = right;
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
  if (left < j) {
    HoareSort(arr, left, j);
  }
  if (i < right) {
    HoareSort(arr, i, right);
  }
}

void ShekhirevHoareBatcherSortSEQ::BatcherMerge(std::vector<int> &arr, int left, int right, int step) {
  int n = right - left + 1;
  if (n <= step) {
    return;
  }
  BatcherMerge(arr, left, right, step * 2);
  BatcherMerge(arr, left + step, right, step * 2);
  for (int i = left + step; i + step <= right; i += step * 2) {
    if (arr[i] > arr[i + step]) {
      std::swap(arr[i], arr[i + step]);
    }
  }
}

bool ShekhirevHoareBatcherSortSEQ::RunImpl() {
  auto &data = GetOutput();
  if (data.size() <= 1) {
    return true;
  }

  int n = static_cast<int>(data.size());
  int mid = n / 2;

  HoareSort(data, 0, mid - 1);
  HoareSort(data, mid, n - 1);

  BatcherMerge(data, 0, n - 1, 1);

  return true;
}

bool ShekhirevHoareBatcherSortSEQ::PostProcessingImpl() {
  size_t original_size = GetInput().size();
  GetOutput().resize(original_size);
  return true;
}

}  // namespace shekhirev_v_hoare_batcher_sort_seq
