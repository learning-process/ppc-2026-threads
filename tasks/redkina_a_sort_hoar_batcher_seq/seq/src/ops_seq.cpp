// ops_seq.cpp
#include "redkina_a_sort_hoar_batcher_seq/seq/include/ops_seq.hpp"

#include <algorithm>
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

// Итеративная версия OddEvenMerge без рекурсии
std::vector<int> OddEvenMerge(const std::vector<int> &left, const std::vector<int> &right) {
  struct Task {
    std::vector<int> left_vec;
    std::vector<int> right_vec;
    int stage;  // 0 = split, 1 = merge even, 2 = merge odd, 3 = combine
    std::vector<int> even_result;
    std::vector<int> odd_result;
  };
  std::stack<Task> stack;
  stack.push({left, right, 0, {}, {}});
  std::vector<int> final_result;

  while (!stack.empty()) {
    Task &task = stack.top();
    if (task.left_vec.empty() && task.right_vec.empty() && task.stage == 3) {
      final_result = Interleave(task.even_result, task.odd_result);
      FinalCompare(final_result);
      stack.pop();
      continue;
    }

    // Базовые случаи
    if (task.stage == 0) {
      if (task.left_vec.empty()) {
        final_result = task.right_vec;
        stack.pop();
        continue;
      }
      if (task.right_vec.empty()) {
        final_result = task.left_vec;
        stack.pop();
        continue;
      }
      if (task.left_vec.size() == 1 && task.right_vec.size() == 1) {
        if (task.left_vec[0] <= task.right_vec[0]) {
          final_result = {task.left_vec[0], task.right_vec[0]};
        } else {
          final_result = {task.right_vec[0], task.left_vec[0]};
        }
        stack.pop();
        continue;
      }

      // Разделяем на чётные/нечётные
      std::vector<int> left_even, left_odd, right_even, right_odd;
      SplitEvenOdd(task.left_vec, left_even, left_odd);
      SplitEvenOdd(task.right_vec, right_even, right_odd);
      task.stage = 1;
      // Кладём подзадачи в стек (обратный порядок для правильного выполнения)
      stack.push({left_odd, right_odd, 0, {}, {}});    // сначала обработаем нечётные
      stack.push({left_even, right_even, 0, {}, {}});  // потом чётные
    } else if (task.stage == 1) {
      // Обработаны чётные (результат в final_result)
      task.even_result = final_result;
      task.stage = 2;
    } else if (task.stage == 2) {
      // Обработаны нечётные (результат в final_result)
      task.odd_result = final_result;
      task.stage = 3;
    }
  }
  return final_result;
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
      int mid = l + (len / 2);
      top.stage = 1;
      stack.push({l, mid, 0});
    } else if (top.stage == 1) {
      int mid = l + (len / 2);
      top.stage = 2;
      stack.push({mid, r, 0});
    } else {
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
