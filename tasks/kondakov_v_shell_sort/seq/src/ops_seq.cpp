#include "kondakov_v_shell_sort/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "kondakov_v_shell_sort/common/include/common.hpp"

namespace kondakov_v_shell_sort {

namespace {

void ShellSort(std::vector<int> &data) {
  const size_t n = data.size();
  for (size_t gap = n / 2; gap > 0; gap /= 2) {
    for (size_t i = gap; i < n; ++i) {
      int value = data[i];
      size_t j = i;
      while (j >= gap && data[j - gap] > value) {
        data[j] = data[j - gap];
        j -= gap;
      }
      data[j] = value;
    }
  }
}

void SimpleMerge(const std::vector<int> &left, const std::vector<int> &right, std::vector<int> &result) {
  size_t i = 0;
  size_t j = 0;
  size_t k = 0;

  while (i < left.size() && j < right.size()) {
    if (left[i] <= right[j]) {
      result[k++] = left[i++];
    } else {
      result[k++] = right[j++];
    }
  }

  while (i < left.size()) {
    result[k++] = left[i++];
  }

  while (j < right.size()) {
    result[k++] = right[j++];
  }
}

}  // namespace

KondakovVShellSortSEQ::KondakovVShellSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool KondakovVShellSortSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool KondakovVShellSortSEQ::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool KondakovVShellSortSEQ::RunImpl() {
  std::vector<int> &data = GetOutput();
  if (data.size() <= 1) {
    return true;
  }

  const size_t mid = data.size() / 2;
  std::vector<int> left_part(data.begin(), data.begin() + static_cast<std::ptrdiff_t>(mid));
  std::vector<int> right_part(data.begin() + static_cast<std::ptrdiff_t>(mid), data.end());

  ShellSort(left_part);
  ShellSort(right_part);
  SimpleMerge(left_part, right_part, data);

  return std::ranges::is_sorted(data);
}

bool KondakovVShellSortSEQ::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace kondakov_v_shell_sort
