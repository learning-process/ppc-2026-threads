#include "trofimov_n_hoar_sort_batcher/seq/include/ops_seq.hpp"

#include <algorithm>
#include <vector>

namespace trofimov_n_hoar_sort_batcher {

namespace {

int HoarePartition(std::vector<int> &arr, int left, int right) {
  int pivot = arr[left + (right - left) / 2];
  int i = left - 1;
  int j = right + 1;

  while (true) {
    do {
      ++i;
    } while (arr[i] < pivot);

    do {
      --j;
    } while (arr[j] > pivot);

    if (i >= j) {
      return j;
    }

    std::swap(arr[i], arr[j]);
  }
}

void HoareQuickSort(std::vector<int> &arr, int left, int right) {
  if (left >= right) {
    return;
  }

  int pivot_index = HoarePartition(arr, left, right);

  HoareQuickSort(arr, left, pivot_index);
  HoareQuickSort(arr, pivot_index + 1, right);
}

void CompareExchange(int &a, int &b) {
  if (a > b) {
    std::swap(a, b);
  }
}

void OddEvenMerge(std::vector<int> &arr, int left, int right, int step) {
  int dist = step * 2;
  if (dist < right - left + 1) {
    OddEvenMerge(arr, left, right, dist);
    OddEvenMerge(arr, left + step, right, dist);

    for (int i = left + step; i + step <= right; i += dist) {
      CompareExchange(arr[i], arr[i + step]);
    }
  } else {
    if (left + step <= right) {
      CompareExchange(arr[left], arr[left + step]);
    }
  }
}

void BatcherMergeSort(std::vector<int> &arr, int left, int right) {
  int size = right - left + 1;
  if (size <= 1) {
    return;
  }

  int mid = left + size / 2 - 1;

  BatcherMergeSort(arr, left, mid);
  BatcherMergeSort(arr, mid + 1, right);

  OddEvenMerge(arr, left, right, 1);
}

}  // namespace

TrofimovNHoarSortBatcherSEQ::TrofimovNHoarSortBatcherSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TrofimovNHoarSortBatcherSEQ::ValidationImpl() {
  return true;
}

bool TrofimovNHoarSortBatcherSEQ::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool TrofimovNHoarSortBatcherSEQ::RunImpl() {
  auto &data = GetOutput();

  if (data.size() > 1) {
    HoareQuickSort(data, 0, static_cast<int>(data.size()) - 1);
    BatcherMergeSort(data, 0, static_cast<int>(data.size()) - 1);
  }

  return true;
}

bool TrofimovNHoarSortBatcherSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace trofimov_n_hoar_sort_batcher
