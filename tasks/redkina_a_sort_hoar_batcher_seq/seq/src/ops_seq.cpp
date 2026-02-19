// ops_seq.cpp
#include "redkina_a_sort_hoar_batcher_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <stack>
#include <utility>
#include <vector>

#include "redkina_a_sort_hoar_batcher_seq/common/include/common.hpp"

namespace redkina_a_sort_hoar_batcher_seq {

namespace {
constexpr int kQuickSortThreshold = 16;

int Partition(std::vector<int> &arr, int left, int right, int pivot) {
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
  return j;
}

void ProcessOneIteration(std::vector<int> &arr, int l, int r, std::stack<std::pair<int, int>> &stack) {
  if (l >= r) {
    return;
  }
  int pivot = arr[(l + r) / 2];
  int j = Partition(arr, l, r, pivot);
  if (l < j) {
    stack.emplace(l, j);
  }
  if (j + 1 < r) {
    stack.emplace(j + 1, r);
  }
}

void QuickSortIterative(std::vector<int> &arr, int left, int right) {
  std::stack<std::pair<int, int>> stack;
  stack.emplace(left, right);
  while (!stack.empty()) {
    auto [l, r] = stack.top();
    stack.pop();
    ProcessOneIteration(arr, l, r, stack);
  }
}

void SplitEvenOdd(const std::vector<int> &src, std::vector<int> &even, std::vector<int> &odd) {
  even.clear();
  odd.clear();
  for (size_t i = 0; i < src.size(); ++i) {
    if (i % 2 == 0) {
      even.push_back(src[i]);
    } else {
      odd.push_back(src[i]);
    }
  }
}

std::vector<int> Interleave(const std::vector<int> &even, const std::vector<int> &odd) {
  std::vector<int> result;
  size_t max_size = std::max(even.size(), odd.size());
  result.reserve(even.size() + odd.size());
  for (size_t i = 0; i < max_size; ++i) {
    if (i < even.size()) {
      result.push_back(even[i]);
    }
    if (i < odd.size()) {
      result.push_back(odd[i]);
    }
  }
  return result;
}

void FinalCompare(std::vector<int> &arr) {
  for (size_t i = 1; i + 1 < arr.size(); i += 2) {
    if (arr[i] > arr[i + 1]) {
      std::swap(arr[i], arr[i + 1]);
    }
  }
}

std::vector<int> HandleBaseCases(const std::vector<int> &left, const std::vector<int> &right) {
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
  return {};
}

std::vector<int> OddEvenMerge(const std::vector<int> &left, const std::vector<int> &right) {
  std::vector<int> base = HandleBaseCases(left, right);
  if (!base.empty()) {
    return base;
  }

  std::vector<int> left_even;
  std::vector<int> left_odd;
  std::vector<int> right_even;
  std::vector<int> right_odd;
  SplitEvenOdd(left, left_even, left_odd);
  SplitEvenOdd(right, right_even, right_odd);

  std::vector<int> merged_even = OddEvenMerge(left_even, right_even);
  std::vector<int> merged_odd = OddEvenMerge(left_odd, right_odd);

  std::vector<int> result = Interleave(merged_even, merged_odd);
  FinalCompare(result);
  return result;
}

struct SortRange {
  int left;
  int right;
  int stage;  // 0 = need split, 1 = need merge after left done, 2 = need merge after both done
};

void HybridSort(std::vector<int> &arr, int left, int right) {
  std::stack<SortRange> stack;
  stack.push({left, right, 0});
  while (!stack.empty()) {
    SortRange &top = stack.top();
    int l = top.left;
    int r = top.right;
    int len = r - l;

    if (len <= 1) {
      stack.pop();
      continue;
    }

    if (len <= kQuickSortThreshold) {
      QuickSortIterative(arr, l, r - 1);
      stack.pop();
      continue;
    }

    if (top.stage == 0) {
      // first visit: split and push left half
      int mid = l + (len / 2);
      top.stage = 1;
      stack.push({l, mid, 0});
    } else if (top.stage == 1) {
      // after left half done, push right half
      int mid = l + (len / 2);
      top.stage = 2;
      stack.push({mid, r, 0});
    } else {
      // both halves done, merge
      int mid = l + (len / 2);
      std::vector<int> left_part(arr.begin() + l, arr.begin() + mid);
      std::vector<int> right_part(arr.begin() + mid, arr.begin() + r);
      std::vector<int> merged = OddEvenMerge(left_part, right_part);
      std::ranges::copy(merged, arr.begin() + l);
      stack.pop();
    }
  }
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
