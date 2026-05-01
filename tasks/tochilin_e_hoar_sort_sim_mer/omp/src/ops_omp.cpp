#include "tochilin_e_hoar_sort_sim_mer/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

#include "tochilin_e_hoar_sort_sim_mer/common/include/common.hpp"

namespace tochilin_e_hoar_sort_sim_mer {

TochilinEHoarSortSimMerOMP::TochilinEHoarSortSimMerOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TochilinEHoarSortSimMerOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool TochilinEHoarSortSimMerOMP::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

std::pair<int, int> TochilinEHoarSortSimMerOMP::Partition(std::vector<int> &arr, int l, int r) {
  int i = l;
  int j = r;
  const int pivot = arr[(l + r) / 2];

  while (i <= j) {
    while (arr[i] < pivot) {
      ++i;
    }
    while (arr[j] > pivot) {
      --j;
    }
    if (i <= j) {
      std::swap(arr[i], arr[j]);
      ++i;
      --j;
    }
  }

  return {i, j};
}

void TochilinEHoarSortSimMerOMP::QuickSortOMP(std::vector<int> &arr, int low, int high, int depth_limit) {
  std::vector<std::pair<int, int>> stack;
  stack.emplace_back(low, high);

  while (!stack.empty()) {
    const auto [l0, r0] = stack.back();
    stack.pop_back();

    int l = l0;
    int r = r0;

    if (l >= r) {
      continue;
    }

    const std::pair<int, int> bounds = Partition(arr, l, r);
    int i = bounds.first;
    int j = bounds.second;

    const bool spawn_tasks = depth_limit > 0;

    if (spawn_tasks) {
#pragma omp task default(none) shared(arr) firstprivate(l, j, depth_limit)
      QuickSortOMP(arr, l, j, depth_limit - 1);

#pragma omp task default(none) shared(arr) firstprivate(i, r, depth_limit)
      QuickSortOMP(arr, i, r, depth_limit - 1);
    } else {
      if (l < j) {
        stack.emplace_back(l, j);
      }
      if (i < r) {
        stack.emplace_back(i, r);
      }
    }
  }

#pragma omp taskwait
}

std::vector<int> TochilinEHoarSortSimMerOMP::MergeSortedVectors(const std::vector<int> &a, const std::vector<int> &b) {
  std::vector<int> result;
  result.reserve(a.size() + b.size());

  std::ranges::merge(a, b, std::back_inserter(result));
  return result;
}

bool TochilinEHoarSortSimMerOMP::RunImpl() {
  auto &data = GetOutput();

  if (data.empty()) {
    return false;
  }

  const auto mid = static_cast<std::vector<int>::difference_type>(data.size() / 2);

  std::vector<int> left(data.begin(), data.begin() + mid);
  std::vector<int> right(data.begin() + mid, data.end());

#pragma omp parallel default(none) shared(left, right)
  {
#pragma omp single
    {
      QuickSortOMP(left, 0, static_cast<int>(left.size()) - 1, 3);
      QuickSortOMP(right, 0, static_cast<int>(right.size()) - 1, 3);
    }
  }

  data = MergeSortedVectors(left, right);
  return true;
}

bool TochilinEHoarSortSimMerOMP::PostProcessingImpl() {
  return std::ranges::is_sorted(GetOutput());
}

}  // namespace tochilin_e_hoar_sort_sim_mer
