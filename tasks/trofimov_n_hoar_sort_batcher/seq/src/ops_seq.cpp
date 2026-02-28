#include "trofimov_n_hoar_sort_batcher/seq/include/ops_seq.hpp"

#include <algorithm>
#include <stack>
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

void CompareExchange(int &a, int &b) {
  if (a > b) {
    std::swap(a, b);
  }
}

void OddEvenMergeIter(std::vector<int> &arr, int left, int right) {
  int n = right - left + 1;
  for (int step = 1; step < n; step *= 2) {
    for (int i = left; i + step <= right; i += step * 2) {
      CompareExchange(arr[i], arr[i + step]);
    }
  }
}

void QuickBatcherIterative(std::vector<int> &arr, int left, int right) {
  std::stack<std::pair<int, int>> stk;
  stk.emplace(left, right);

  while (!stk.empty()) {
    auto [l, r] = stk.top();
    stk.pop();
    if (l >= r) {
      continue;
    }

    int pivot_index = HoarePartition(arr, l, r);

    if (pivot_index - l > r - pivot_index - 1) {
      stk.emplace(l, pivot_index);
      stk.emplace(pivot_index + 1, r);
    } else {
      stk.emplace(pivot_index + 1, r);
      stk.emplace(l, pivot_index);
    }

    OddEvenMergeIter(arr, l, r);
  }
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
    QuickBatcherIterative(data, 0, static_cast<int>(data.size()) - 1);
  }
  return true;
}

bool TrofimovNHoarSortBatcherSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace trofimov_n_hoar_sort_batcher
