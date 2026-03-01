#include "shekhirev_v_hoare_batcher_sort_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <climits>
#include <cstddef>
#include <tuple>
#include <utility>
#include <vector>

#include "shekhirev_v_hoare_batcher_sort_seq/common/include/common.hpp"

namespace shekhirev_v_hoare_batcher_sort_seq {

namespace {

int HoarePartition(std::vector<int> &arr, int left, int right) {
  int pivot = arr[left + ((right - left) / 2)];
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
  return i;
}

void BatcherComp(std::vector<int> &arr, int left, int right, int r, int m) {
  for (int i = left + r; i + r <= right; i += m) {
    if (arr[i] > arr[i + r]) {
      std::swap(arr[i], arr[i + r]);
    }
  }
}

}  // namespace

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
  if (in.empty()) {
    GetOutput().clear();
    return true;
  }

  size_t original_size = in.size();
  size_t p2 = 1;
  while (p2 < original_size) {
    p2 *= 2;
  }

  GetOutput() = in;
  GetOutput().resize(p2, INT_MAX);
  return true;
}

void ShekhirevHoareBatcherSortSEQ::HoareSort(std::vector<int> &arr, int start_left, int start_right) {
  std::vector<std::pair<int, int>> stk;
  stk.reserve(32);
  stk.emplace_back(start_left, start_right);

  while (!stk.empty()) {
    auto [left, right] = stk.back();
    stk.pop_back();

    if (left >= right) {
      continue;
    }

    int partition_index = HoarePartition(arr, left, right);
    if (partition_index < right) {
      stk.emplace_back(partition_index, right);
    }
    if (left < partition_index - 1) {
      stk.emplace_back(left, partition_index - 1);
    }
  }
}

void ShekhirevHoareBatcherSortSEQ::BatcherMerge(std::vector<int> &arr, int start_left, int start_right, int start_r) {
  std::vector<std::tuple<int, int, int, int>> stk;
  stk.reserve(32);
  stk.emplace_back(start_left, start_right, start_r, 0);

  while (!stk.empty()) {
    auto [left, right, r, stage] = stk.back();
    stk.pop_back();

    int n = right - left + 1;
    int m = r * 2;

    if (m < n) {
      if (stage == 0) {
        stk.emplace_back(left, right, r, 1);
        stk.emplace_back(left + r, right, m, 0);
        stk.emplace_back(left, right, m, 0);
      } else {
        BatcherComp(arr, left, right, r, m);
      }
    } else {
      if ((left + r) <= right) {
        if (arr[left] > arr[left + r]) {
          std::swap(arr[left], arr[left + r]);
        }
      }
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
