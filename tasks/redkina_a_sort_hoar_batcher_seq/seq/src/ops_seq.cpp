// ops_seq.cpp
#include "redkina_a_sort_hoar_batcher_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

namespace redkina_a_sort_hoar_batcher_seq {

namespace {
constexpr int kQuickSortThreshold = 16;

void QuickSort(std::vector<int> &arr, int left, int right) {
  if (left >= right) {
    return;
  }
  int pivot = arr[(left + right) / 2];
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
  if (left < j) {
    QuickSort(arr, left, j);
  }
  if (i < right) {
    QuickSort(arr, i, right);
  }
}

std::vector<int> OddEvenMerge(const std::vector<int> &left, const std::vector<int> &right) {
  if (left.empty()) {
    return right;
  }
  if (right.empty()) {
    return left;
  }
  if (left.size() == 1 && right.size() == 1) {
    if (left[0] <= right[0]) {
      return {left[0], right[0]};
    }
    return {right[0], left[0]};
  }

  std::vector<int> left_even;
  std::vector<int> left_odd;
  std::vector<int> right_even;
  std::vector<int> right_odd;
  for (size_t i = 0; i < left.size(); ++i) {
    if (i % 2 == 0) {
      left_even.push_back(left[i]);
    } else {
      left_odd.push_back(left[i]);
    }
  }
  for (size_t i = 0; i < right.size(); ++i) {
    if (i % 2 == 0) {
      right_even.push_back(right[i]);
    } else {
      right_odd.push_back(right[i]);
    }
  }

  std::vector<int> merged_even = OddEvenMerge(left_even, right_even);
  std::vector<int> merged_odd = OddEvenMerge(left_odd, right_odd);

  std::vector<int> result;
  size_t max_size = std::max(merged_even.size(), merged_odd.size());
  for (size_t i = 0; i < max_size; ++i) {
    if (i < merged_even.size()) {
      result.push_back(merged_even[i]);
    }
    if (i < merged_odd.size()) {
      result.push_back(merged_odd[i]);
    }
  }

  for (size_t i = 1; i + 1 < result.size(); i += 2) {
    if (result[i] > result[i + 1]) {
      std::swap(result[i], result[i + 1]);
    }
  }
  return result;
}

void HybridSort(std::vector<int> &arr, int left, int right) {
  int len = right - left;
  if (len <= 1) {
    return;
  }
  if (len <= kQuickSortThreshold) {
    QuickSort(arr, left, right - 1);
    return;
  }
  int mid = left + len / 2;
  HybridSort(arr, left, mid);
  HybridSort(arr, mid, right);

  std::vector<int> left_part(arr.begin() + left, arr.begin() + mid);
  std::vector<int> right_part(arr.begin() + mid, arr.begin() + right);

  std::vector<int> merged = OddEvenMerge(left_part, right_part);

  std::ranges::copy(merged, arr.begin() + left);
}

}  // namespace

RedkinaASortHoarBatcherSEQ::RedkinaASortHoarBatcherSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool RedkinaASortHoarBatcherSEQ::ValidationImpl() {
  return true;
}

bool RedkinaASortHoarBatcherSEQ::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool RedkinaASortHoarBatcherSEQ::RunImpl() {
  auto &out = GetOutput();
  if (out.empty()) {
    return true;
  }
  HybridSort(out, 0, static_cast<int>(out.size()));
  return true;
}

bool RedkinaASortHoarBatcherSEQ::PostProcessingImpl() {
  return GetOutput().size() == GetInput().size();
}

}  // namespace redkina_a_sort_hoar_batcher_seq
