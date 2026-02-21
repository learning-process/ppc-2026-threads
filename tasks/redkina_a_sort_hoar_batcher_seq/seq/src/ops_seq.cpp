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
  const size_t max_size = std::max(even.size(), odd.size());
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

inline std::vector<int> OddEvenMerge(const std::vector<int> &left, const std::vector<int> &right) {
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

void HybridSort(std::vector<int> &arr, const int left, const int right) {
  std::stack<SortRange> stack;
  stack.push({left, right, 0});
  while (!stack.empty()) {
    SortRange &top = stack.top();
    const int l = top.left;
    const int r = top.right;
    const int len = r - l;

    if (len <= 1) {
      stack.pop();
      continue;
    }

    if (len <= kQuickSortThreshold) {
      std::sort(arr.begin() + l, arr.begin() + r);
      stack.pop();
      continue;
    }

    if (top.stage == 0) {
      const int mid = l + (len / 2);
      top.stage = 1;
      stack.push({l, mid, 0});
    } else if (top.stage == 1) {
      const int mid = l + (len / 2);
      top.stage = 2;
      stack.push({mid, r, 0});
    } else {
      const int mid = l + (len / 2);
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
