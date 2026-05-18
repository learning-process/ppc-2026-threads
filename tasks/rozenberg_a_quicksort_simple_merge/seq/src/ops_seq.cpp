#include "rozenberg_a_quicksort_simple_merge/seq/include/ops_seq.hpp"

#include <stack>
#include <utility>
#include <vector>

#include "rozenberg_a_quicksort_simple_merge/common/include/common.hpp"

namespace rozenberg_a_quicksort_simple_merge {

RozenbergAQuicksortSimpleMergeSEQ::RozenbergAQuicksortSimpleMergeSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());

  InType empty;
  GetInput().swap(empty);

  for (const auto &elem : in) {
    GetInput().push_back(elem);
  }

  GetOutput().clear();
}

bool RozenbergAQuicksortSimpleMergeSEQ::ValidationImpl() {
  return (!(GetInput().empty())) && (GetOutput().empty());
}

bool RozenbergAQuicksortSimpleMergeSEQ::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return GetOutput().size() == GetInput().size();
}

std::pair<int, int> RozenbergAQuicksortSimpleMergeSEQ::Partition(InType &data, int left, int right) {
  const int pivot = data[left + ((right - left) / 2)];
  int i = left;
  int j = right;

  while (i <= j) {
    while (data[i] < pivot) {
      i++;
    }
    while (data[j] > pivot) {
      j--;
    }

    if (i <= j) {
      std::swap(data[i], data[j]);
      i++;
      j--;
    }
  }
  return {i, j};
}

void RozenbergAQuicksortSimpleMergeSEQ::PushSubarrays(std::stack<std::pair<int, int>> &stack, int left, int right,
                                                      int i, int j) {
  if (j - left > right - i) {
    if (left < j) {
      stack.emplace(left, j);
    }
    if (i < right) {
      stack.emplace(i, right);
    }
  } else {
    if (i < right) {
      stack.emplace(i, right);
    }
    if (left < j) {
      stack.emplace(left, j);
    }
  }
}

void RozenbergAQuicksortSimpleMergeSEQ::Quicksort(InType &data) {
  if (data.size() < 2) {
    return;
  }

  std::stack<std::pair<int, int>> stack;

  stack.emplace(0, data.size() - 1);

  while (!stack.empty()) {
    const auto [left, right] = stack.top();
    stack.pop();

    if (left < right) {
      const auto [i, j] = Partition(data, left, right);
      PushSubarrays(stack, left, right, i, j);
    }
  }
}

bool RozenbergAQuicksortSimpleMergeSEQ::RunImpl() {
  InType data = GetInput();
  Quicksort(data);
  GetOutput() = data;
  return true;
}

bool RozenbergAQuicksortSimpleMergeSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace rozenberg_a_quicksort_simple_merge
