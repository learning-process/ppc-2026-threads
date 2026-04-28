#include "trofimov_n_hoar_sort_batcher/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <vector>

#include "oneapi/tbb/task_group.h"
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

void QuickSortTbbTask(std::vector<int> &arr, int left, int right, int depth_limit) {
  if (left >= right) {
    return;
  }

  constexpr int kSequentialThreshold = 1024;

  if ((right - left) < kSequentialThreshold || depth_limit <= 0) {
    std::sort(arr.begin() + left, arr.begin() + right + 1);
    return;
  }

  const int split = HoarePartition(arr, left, right);

  oneapi::tbb::task_group tg;
  if ((split - left) > kSequentialThreshold) {
    tg.run([&arr, left, split, depth_limit]() { QuickSortTbbTask(arr, left, split, depth_limit - 1); });
  } else {
    std::sort(arr.begin() + left, arr.begin() + split + 1);
  }

  if ((right - (split + 1)) > kSequentialThreshold) {
    tg.run([&arr, split, right, depth_limit]() { QuickSortTbbTask(arr, split + 1, right, depth_limit - 1); });
  } else {
    std::sort(arr.begin() + split + 1, arr.begin() + right + 1);
  }

  tg.wait();
}

}  // namespace

TrofimovNHoarSortBatcherTBB::TrofimovNHoarSortBatcherTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TrofimovNHoarSortBatcherTBB::ValidationImpl() {
  return true;
}

bool TrofimovNHoarSortBatcherTBB::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool TrofimovNHoarSortBatcherTBB::RunImpl() {
  auto &data = GetOutput();

  if (data.size() > 1) {
    QuickSortTbbTask(data, 0, static_cast<int>(data.size()) - 1, 4);
  }

  return true;
}

bool TrofimovNHoarSortBatcherTBB::PostProcessingImpl() {
  return true;
}

}  // namespace trofimov_n_hoar_sort_batcher
