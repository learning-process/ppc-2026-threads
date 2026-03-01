#include "likhanov_m_hoare_batcher_sort/seq/include/ops_seq.hpp"

#include <algorithm>
#include <stack>
#include <utility>
#include <vector>

#include "likhanov_m_hoare_batcher_sort/common/include/common.hpp"

namespace likhanov_m_hoare_batcher_sort {

namespace {

int Partition(std::vector<int> &arr, int left, int right) {
  const int pivot = arr[(left + right) / 2];
  int i = left;
  int j = right;

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

  return i;
}

void HoareSort(std::vector<int> &arr, int left, int right) {
  if (left >= right) {
    return;
  }

  std::stack<std::pair<int, int>> stack;
  stack.emplace(left, right);

  while (!stack.empty()) {
    const auto [l, r] = stack.top();
    stack.pop();

    if (l >= r) {
      continue;
    }

    const int index = Partition(arr, l, r);

    if (l < index - 1) {
      stack.emplace(l, index - 1);
    }
    if (index < r) {
      stack.emplace(index, r);
    }
  }
}

void OddEvenMerge(std::vector<int> &arr, int left, int mid, int right) {
  const int n = right - left + 1;
  if (n <= 1) {
    return;
  }

  for (int step = mid - left; step > 0; step /= 2) {
    for (int i = left; i + step <= right; ++i) {
      if (((i - left) / step) % 2 == 0) {
        if (arr[i] > arr[i + step]) {
          std::swap(arr[i], arr[i + step]);
        }
      }
    }
  }
}

void BatcherMerge(std::vector<int> &arr, int left, int mid, int right) {
  std::inplace_merge(arr.begin() + left, arr.begin() + mid + 1, arr.begin() + right + 1);

  OddEvenMerge(arr, left, mid, right);
}

}  // namespace

LikhanovMHoareBatcherSortSEQ::LikhanovMHoareBatcherSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool LikhanovMHoareBatcherSortSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool LikhanovMHoareBatcherSortSEQ::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool LikhanovMHoareBatcherSortSEQ::RunImpl() {
  auto &data = GetOutput();
  const int n = static_cast<int>(data.size());

  if (n <= 1) {
    return true;
  }

  const int mid = (n / 2) - 1;

  HoareSort(data, 0, mid);
  HoareSort(data, mid + 1, n - 1);

  BatcherMerge(data, 0, mid, n - 1);

  return true;
}

bool LikhanovMHoareBatcherSortSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace likhanov_m_hoare_batcher_sort
