#include "trofimov_n_hoar_sort_batcher/stl/include/ops_stl.hpp"

#include <algorithm>
#include <thread>
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

void QuickSortStlTask(std::vector<int> &arr, int left, int right, int depth_limit) {
  if (left >= right) {
    return;
  }

  constexpr int kSequentialThreshold = 1024;

  if ((right - left) < kSequentialThreshold || depth_limit <= 0) {
    std::sort(arr.begin() + left, arr.begin() + right + 1);
    return;
  }

  const int split = HoarePartition(arr, left, right);

  std::thread left_thread([&arr, left, split, depth_limit]() { QuickSortStlTask(arr, left, split, depth_limit - 1); });
  QuickSortStlTask(arr, split + 1, right, depth_limit - 1);
  left_thread.join();
}

}  // namespace

TrofimovNHoarSortBatcherSTL::TrofimovNHoarSortBatcherSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TrofimovNHoarSortBatcherSTL::ValidationImpl() {
  return true;
}

bool TrofimovNHoarSortBatcherSTL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool TrofimovNHoarSortBatcherSTL::RunImpl() {
  auto &data = GetOutput();

  if (data.size() > 1) {
    QuickSortStlTask(data, 0, static_cast<int>(data.size()) - 1, 4);
  }

  return true;
}

bool TrofimovNHoarSortBatcherSTL::PostProcessingImpl() {
  return true;
}

}  // namespace trofimov_n_hoar_sort_batcher
