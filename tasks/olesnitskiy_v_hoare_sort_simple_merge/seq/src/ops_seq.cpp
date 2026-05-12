#include "olesnitskiy_v_hoare_sort_simple_merge/seq/include/ops_seq.hpp"

#include <algorithm>
#include <stack>
#include <utility>
#include <vector>

#include "olesnitskiy_v_hoare_sort_simple_merge/common/include/common.hpp"

namespace olesnitskiy_v_hoare_sort_simple_merge {

OlesnitskiyVHoareSortSimpleMergeSEQ::OlesnitskiyVHoareSortSimpleMergeSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

int OlesnitskiyVHoareSortSimpleMergeSEQ::HoarePartition(std::vector<int> &values, int left, int right) {
  const int pivot = values[left + ((right - left) / 2)];
  int i = left - 1;
  int j = right + 1;

  while (true) {
    ++i;
    while (values[i] < pivot) {
      ++i;
    }

    --j;
    while (values[j] > pivot) {
      --j;
    }

    if (i >= j) {
      return j;
    }

    std::swap(values[i], values[j]);
  }
}

void OlesnitskiyVHoareSortSimpleMergeSEQ::HoareQuickSort(std::vector<int> &values, int left, int right) {
  std::stack<std::pair<int, int>> ranges;
  ranges.emplace(left, right);

  while (!ranges.empty()) {
    auto [current_left, current_right] = ranges.top();
    ranges.pop();

    if (current_left >= current_right) {
      continue;
    }

    const int partition_index = HoarePartition(values, current_left, current_right);

    if ((partition_index - current_left) > (current_right - (partition_index + 1))) {
      ranges.emplace(current_left, partition_index);
      ranges.emplace(partition_index + 1, current_right);
    } else {
      ranges.emplace(partition_index + 1, current_right);
      ranges.emplace(current_left, partition_index);
    }
  }
}

bool OlesnitskiyVHoareSortSimpleMergeSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool OlesnitskiyVHoareSortSimpleMergeSEQ::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool OlesnitskiyVHoareSortSimpleMergeSEQ::RunImpl() {
  std::vector<int> &values = GetOutput();
  const int n = static_cast<int>(values.size());
  if (n <= 1) {
    return true;
  }

  HoareQuickSort(values, 0, n - 1);
  return std::ranges::is_sorted(values);
}

bool OlesnitskiyVHoareSortSimpleMergeSEQ::PostProcessingImpl() {
  return !GetOutput().empty() && std::ranges::is_sorted(GetOutput());
}

}  // namespace olesnitskiy_v_hoare_sort_simple_merge
