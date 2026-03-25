#include "olesnitskiy_v_hoare_sort_simple_merge_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <stack>
#include <utility>
#include <vector>

#include "olesnitskiy_v_hoare_sort_simple_merge_seq/common/include/common.hpp"

namespace olesnitskiy_v_hoare_sort_simple_merge_seq {

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

void OlesnitskiyVHoareSortSimpleMergeSEQ::Merge(std::vector<int> &values, int left, int mid, int right) {
  std::vector<int> merged;
  const int merged_size = (right - left) + 1;
  merged.reserve(static_cast<std::size_t>(merged_size));

  int left_index = left;
  int right_index = mid + 1;

  while (left_index <= mid && right_index <= right) {
    if (values[left_index] <= values[right_index]) {
      merged.push_back(values[left_index++]);
    } else {
      merged.push_back(values[right_index++]);
    }
  }

  while (left_index <= mid) {
    merged.push_back(values[left_index++]);
  }

  while (right_index <= right) {
    merged.push_back(values[right_index++]);
  }

  for (std::size_t idx = 0; idx < merged.size(); ++idx) {
    values[static_cast<std::size_t>(left) + idx] = merged[idx];
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

}  // namespace olesnitskiy_v_hoare_sort_simple_merge_seq
