#include "trofimov_n_hoar_sort_batcher/seq/include/ops_seq.hpp"

#include <algorithm>
#include <vector>

#include "trofimov_n_hoar_sort_batcher/common/include/common.hpp"

namespace trofimov_n_hoar_sort_batcher {

namespace {

int HoarePartition(std::vector<int> &arr, int left, int right) {
  int pivot = arr[left + ((right - left) / 2)];
  int i = left - 1;
  int j = right + 1;

  while (true) {
    ++i;
    while (arr[i] < pivot) {
      ++i;
    }

    --j;
    while (arr[j] > pivot) {
      --j;
    }

    if (i >= j) {
      return j;
    }
    std::swap(arr[i], arr[j]);
  }
}

void OddEvenMergeIter(std::vector<int> &arr, int left, int right) {
  int n = right - left + 1;
  for (int step = 1; step < n; step *= 2) {
    for (int i = left; i + step < left + n; i += step * 2) {
      int j = i + step;
      if (j <= right) {
        if (arr[i] > arr[j]) {
          std::swap(arr[i], arr[j]);
        }
      }
    }
  }
}

void QuickBatcherHybrid(std::vector<int> &arr, int left, int right) {
  if (left >= right) {
    return;
  }

  int pivot_index = HoarePartition(arr, left, right);

  QuickBatcherHybrid(arr, left, pivot_index);
  QuickBatcherHybrid(arr, pivot_index + 1, right);

  OddEvenMergeIter(arr, left, right);
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
    QuickBatcherHybrid(data, 0, static_cast<int>(data.size()) - 1);
  }
  return true;
}

bool TrofimovNHoarSortBatcherSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace trofimov_n_hoar_sort_batcher
