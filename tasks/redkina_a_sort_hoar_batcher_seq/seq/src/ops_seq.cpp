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

// Функция обычного слияния двух отсортированных векторов (итеративная, без рекурсии)
std::vector<int> MergeSorted(const std::vector<int> &left, const std::vector<int> &right) {
  std::vector<int> result;
  result.reserve(left.size() + right.size());
  size_t i = 0, j = 0;
  while (i < left.size() && j < right.size()) {
    if (left[i] <= right[j]) {
      result.push_back(left[i]);
      ++i;
    } else {
      result.push_back(right[j]);
      ++j;
    }
  }
  while (i < left.size()) {
    result.push_back(left[i]);
    ++i;
  }
  while (j < right.size()) {
    result.push_back(right[j]);
    ++j;
  }
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
      // Используем обычное слияние вместо рекурсивной сети Бэтчера
      std::vector<int> merged = MergeSorted(left_part, right_part);
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
