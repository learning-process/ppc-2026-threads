#include "trofimov_n_hoar_sort_batcher/omp/include/ops_omp.hpp"

#include <algorithm>
#include <vector>

#include "trofimov_n_hoar_sort_batcher/common/include/common.hpp"

namespace trofimov_n_hoar_sort_batcher {

namespace {

int HoarePartition(std::vector<int> &arr, int left, int right) {
  const int pivot = arr[left + ((right - left) / 2)];
  int i = left - 1;
  int j = right + 1;

  while (true) {
    while (arr[++i] < pivot) {
    }

    while (arr[--j] > pivot) {
    }

    if (i >= j) {
      return j;
    }

    std::swap(arr[i], arr[j]);
  }
}

void QuickSortOmpTask(std::vector<int> &arr, int left, int right, int depth_limit) {
  if (left >= right) {
    return;
  }

  constexpr int kSequentialThreshold = 1024;

  if ((right - left) < kSequentialThreshold || depth_limit <= 0) {
    std::sort(arr.begin() + left, arr.begin() + right + 1);
    return;
  }

  const int split = HoarePartition(arr, left, right);

#pragma omp task default(none) shared(arr) \
    firstprivate(left, split, depth_limit) if ((split - left) > kSequentialThreshold)
  QuickSortOmpTask(arr, left, split, depth_limit - 1);

#pragma omp task default(none) shared(arr) \
    firstprivate(split, right, depth_limit) if ((right - (split + 1)) > kSequentialThreshold)
  QuickSortOmpTask(arr, split + 1, right, depth_limit - 1);

#pragma omp taskwait
}

}  // namespace

TrofimovNHoarSortBatcherOMP::TrofimovNHoarSortBatcherOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TrofimovNHoarSortBatcherOMP::ValidationImpl() {
  return true;
}

bool TrofimovNHoarSortBatcherOMP::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool TrofimovNHoarSortBatcherOMP::RunImpl() {
  auto &data = GetOutput();

  if (data.size() <= 1) {
    return true;
  }

#pragma omp parallel default(none) shared(data)
  {
#pragma omp single nowait
    {
      QuickSortOmpTask(data, 0, static_cast<int>(data.size()) - 1, 4);
    }
  }

  return true;
}

bool TrofimovNHoarSortBatcherOMP::PostProcessingImpl() {
  return true;
}

}  // namespace trofimov_n_hoar_sort_batcher
