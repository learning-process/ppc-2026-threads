#include "yushkova_p_hoare_sorting_simple_merging/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <stack>
#include <utility>
#include <vector>

#include "yushkova_p_hoare_sorting_simple_merging/common/include/common.hpp"

namespace yushkova_p_hoare_sorting_simple_merging {

YushkovaPHoareSortingSimpleMergingSEQ::YushkovaPHoareSortingSimpleMergingSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

int YushkovaPHoareSortingSimpleMergingSEQ::HoarePartition(std::vector<int> &values, int left, int right) {
  int pivot = values[left + ((right - left) / 2)];
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

void YushkovaPHoareSortingSimpleMergingSEQ::HoareQuickSort(std::vector<int> &values, int left, int right) {
  std::stack<std::pair<int, int>> ranges;
  ranges.emplace(left, right);

  while (!ranges.empty()) {
    auto [current_left, current_right] = ranges.top();
    ranges.pop();

    if (current_left >= current_right) {
      continue;
    }

    int partition_index = HoarePartition(values, current_left, current_right);

    if ((partition_index - current_left) > (current_right - (partition_index + 1))) {
      ranges.emplace(current_left, partition_index);
      ranges.emplace(partition_index + 1, current_right);
    } else {
      ranges.emplace(partition_index + 1, current_right);
      ranges.emplace(current_left, partition_index);
    }
  }
}

std::vector<int> YushkovaPHoareSortingSimpleMergingSEQ::SimpleMerge(const std::vector<int> &left,
                                                                    const std::vector<int> &right) {
  std::vector<int> merged;
  merged.reserve(left.size() + right.size());

  size_t left_index = 0;
  size_t right_index = 0;

  while (left_index < left.size() && right_index < right.size()) {
    if (left[left_index] <= right[right_index]) {
      merged.push_back(left[left_index++]);
    } else {
      merged.push_back(right[right_index++]);
    }
  }

  while (left_index < left.size()) {
    merged.push_back(left[left_index++]);
  }

  while (right_index < right.size()) {
    merged.push_back(right[right_index++]);
  }

  return merged;
}

bool YushkovaPHoareSortingSimpleMergingSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool YushkovaPHoareSortingSimpleMergingSEQ::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool YushkovaPHoareSortingSimpleMergingSEQ::RunImpl() {
  std::vector<int> &values = GetOutput();
  if (values.size() <= 1) {
    return true;
  }

  size_t middle = values.size() / 2;
  std::vector<int> left(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(middle));
  std::vector<int> right(values.begin() + static_cast<std::ptrdiff_t>(middle), values.end());

  if (left.size() > 1) {
    HoareQuickSort(left, 0, static_cast<int>(left.size()) - 1);
  }

  if (right.size() > 1) {
    HoareQuickSort(right, 0, static_cast<int>(right.size()) - 1);
  }

  values = SimpleMerge(left, right);
  return std::ranges::is_sorted(values);
}

bool YushkovaPHoareSortingSimpleMergingSEQ::PostProcessingImpl() {
  return !GetOutput().empty() && std::ranges::is_sorted(GetOutput());
}

}  // namespace yushkova_p_hoare_sorting_simple_merging
