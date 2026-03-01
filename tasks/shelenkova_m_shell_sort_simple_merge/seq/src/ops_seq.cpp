#include "shelenkova_m_shell_sort_simple_merge/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "shelenkova_m_shell_sort_simple_merge/common/include/common.hpp"

namespace shelenkova_m_shell_sort_simple_merge {

namespace {

void ShellSort(std::vector<int> &data) {
  const size_t n = data.size();
  if (n <= 1) {
    return;
  }

  for (size_t gap = n / 2; gap > 0; gap /= 2) {
    for (size_t i = gap; i < n; ++i) {
      int temp = data[i];
      size_t j = i;
      while (j >= gap && data[j - gap] > temp) {
        data[j] = data[j - gap];
        j -= gap;
      }
      data[j] = temp;
    }
  }
}

void SimpleMerge(const std::vector<int> &left, const std::vector<int> &right, std::vector<int> &result) {
  size_t left_idx = 0;
  size_t right_idx = 0;
  size_t res_idx = 0;

  while (left_idx < left.size() && right_idx < right.size()) {
    if (left[left_idx] <= right[right_idx]) {
      result[res_idx++] = left[left_idx++];
    } else {
      result[res_idx++] = right[right_idx++];
    }
  }

  while (left_idx < left.size()) {
    result[res_idx++] = left[left_idx++];
  }

  while (right_idx < right.size()) {
    result[res_idx++] = right[right_idx++];
  }
}

}  // namespace

ShelenkovaMShellSortSimpleMergeSEQ::ShelenkovaMShellSortSimpleMergeSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool ShelenkovaMShellSortSimpleMergeSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool ShelenkovaMShellSortSimpleMergeSEQ::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool ShelenkovaMShellSortSimpleMergeSEQ::RunImpl() {
  std::vector<int> &data = GetOutput();
  if (data.size() <= 1) {
    return true;
  }

  const size_t mid = data.size() / 2;
  std::vector<int> left(data.begin(), data.begin() + static_cast<std::ptrdiff_t>(mid));
  std::vector<int> right(data.begin() + static_cast<std::ptrdiff_t>(mid), data.end());

  ShellSort(left);
  ShellSort(right);
  SimpleMerge(left, right, data);

  return std::ranges::is_sorted(data);
}

bool ShelenkovaMShellSortSimpleMergeSEQ::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace shelenkova_m_shell_sort_simple_merge
