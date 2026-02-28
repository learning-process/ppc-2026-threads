#include "tochilin_e_hoar_sort_sim_mer_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

#include "tochilin_e_hoar_sort_sim_mer_seq/common/include/common.hpp"

namespace tochilin_e_hoar_sort_sim_mer_seq {

TochilinEHoarSortSimMerSEQ::TochilinEHoarSortSimMerSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TochilinEHoarSortSimMerSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool TochilinEHoarSortSimMerSEQ::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

std::pair<int, int> TochilinEHoarSortSimMerSEQ::Partition(std::vector<int> &arr,
                                                          int l, int r) {
  int i = l;
  int j = r;
  int pivot = arr[(l + r) / 2];

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

void TochilinEHoarSortSimMerSEQ::QuickSort(std::vector<int> &arr, int low,
                                           int high) {
  if (low >= high) {
    return;
  }

  std::vector<std::pair<int, int>> stack;
  stack.emplace_back(low, high);

  while (!stack.empty()) {
    auto [l, r] = stack.back();
    stack.pop_back();

    auto [i, j] = Partition(arr, l, r);

    if (l < j) {
      stack.emplace_back(l, j);
    }
    if (i < r) {
      stack.emplace_back(i, r);
    }
  }
}

std::vector<int>
TochilinEHoarSortSimMerSEQ::MergeSortedVectors(const std::vector<int> &a,
                                               const std::vector<int> &b) {
  std::vector<int> result;
  result.reserve(a.size() + b.size());
  std::ranges::merge(a, b, std::back_inserter(result));
  return result;
}

bool TochilinEHoarSortSimMerSEQ::RunImpl() {
  auto &data = GetOutput();
  if (data.empty()) {
    return false;
  }

  const auto mid =
      static_cast<std::vector<int>::difference_type>(data.size() / 2);

  std::vector<int> left(data.begin(), data.begin() + mid);
  std::vector<int> right(data.begin() + mid, data.end());

  QuickSort(left, 0, static_cast<int>(left.size()) - 1);
  QuickSort(right, 0, static_cast<int>(right.size()) - 1);

  data = MergeSortedVectors(left, right);

  return true;
}

bool TochilinEHoarSortSimMerSEQ::PostProcessingImpl() {
  return std::is_sorted(GetOutput().begin(), GetOutput().end());
}

} // namespace tochilin_e_hoar_sort_sim_mer_seq
