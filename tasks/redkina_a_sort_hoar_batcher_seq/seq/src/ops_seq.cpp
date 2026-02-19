// ops_seq.cpp
#include "redkina_a_sort_hoar_batcher_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <stack>
#include <vector>

namespace redkina_a_sort_hoar_batcher_seq {

namespace {
constexpr int kQuickSortThreshold = 16;

// Итеративная быстрая сортировка (без рекурсии)
void QuickSortIterative(std::vector<int> &arr, int left, int right) {
  std::stack<std::pair<int, int>> stack;
  stack.push({left, right});
  while (!stack.empty()) {
    auto [l, r] = stack.top();
    stack.pop();
    if (l >= r) {
      continue;
    }
    int pivot = arr[(l + r) / 2];
    int i = l;
    int j = r;
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
    if (l < j) {
      stack.push({l, j});
    }
    if (i < r) {
      stack.push({i, r});
    }
  }
}

// Разделение вектора на элементы с чётными и нечётными индексами
static void SplitEvenOdd(const std::vector<int> &src, std::vector<int> &even, std::vector<int> &odd) {
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

// Чередование двух векторов (сначала чётные, потом нечётные)
static std::vector<int> Interleave(const std::vector<int> &even, const std::vector<int> &odd) {
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

// Финальная фаза слияния Бэтчера: сравнение соседних элементов (начиная с индекса 1)
static void FinalCompare(std::vector<int> &arr) {
  for (size_t i = 1; i + 1 < arr.size(); i += 2) {
    if (arr[i] > arr[i + 1]) {
      std::swap(arr[i], arr[i + 1]);
    }
  }
}

// Рекурсивное слияние Бэтчера (оставлено рекурсивным по алгоритму)
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

  std::vector<int> left_even, left_odd;
  std::vector<int> right_even, right_odd;
  SplitEvenOdd(left, left_even, left_odd);
  SplitEvenOdd(right, right_even, right_odd);

  std::vector<int> merged_even = OddEvenMerge(left_even, right_even);
  std::vector<int> merged_odd = OddEvenMerge(left_odd, right_odd);

  std::vector<int> result = Interleave(merged_even, merged_odd);
  FinalCompare(result);
  return result;
}

// Гибридная сортировка (рекурсивная)
void HybridSort(std::vector<int> &arr, int left, int right) {
  int len = right - left;
  if (len <= 1) {
    return;
  }
  if (len <= kQuickSortThreshold) {
    QuickSortIterative(arr, left, right - 1);
    return;
  }
  int mid = left + (len / 2);
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
