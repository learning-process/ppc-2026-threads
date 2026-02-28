#include "chetverikova_e_shell_sort_simple_merge/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>

#include "chetverikova_e_shell_sort_simple_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace chetverikova_e_shell_sort_simple_merge {

ChetverikovaEShellSortSimpleMergeSEQ::ChetverikovaEShellSortSimpleMergeSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ChetverikovaEShellSortSimpleMergeSEQ::ValidationImpl() {
  return !(GetInput().empty());
}

bool ChetverikovaEShellSortSimpleMergeSEQ::PreProcessingImpl() {
  return true;
}

void ChetverikovaEShellSortSimpleMergeSEQ::ShellSort(std::vector<int> &data) {
  if (data.empty()) {
    return;
  }

  size_t n = data.size();
  for (size_t gap = n / 2; gap > 0; gap /= 2) {
    for (size_t i = gap; i < n; i++) {
      int temp = data[i];
      size_t j;

      for (j = i; j >= gap && data[j - gap] > temp; j -= gap) {
        data[j] = data[j - gap];
      }

      data[j] = temp;
    }
  }
}

std::vector<int> ChetverikovaEShellSortSimpleMergeSEQ::MergeSort(const std::vector<int> &left,
                                                                 const std::vector<int> &right) {
  std::vector<int> result;
  result.reserve(left.size() + right.size());
  size_t i = 0, j = 0;
  while (i < left.size() && j < right.size()) {
    if (left[i] <= right[j]) {
      result.push_back(left[i++]);
    } else {
      result.push_back(right[j++]);
    }
  }
  while (i < left.size()) {
    result.push_back(left[i++]);
  }
  while (j < right.size()) {
    result.push_back(right[j++]);
  }

  return result;
}

bool ChetverikovaEShellSortSimpleMergeSEQ::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  if (input.empty()) {
    GetOutput().clear();
    return true;
  }

  const size_t mid = input.size() / 2;
  std::vector<int> left(input.begin(), input.begin() + static_cast<std::ptrdiff_t>(mid));
  std::vector<int> right(input.begin() + static_cast<std::ptrdiff_t>(mid), input.end());
  ShellSort(left);
  ShellSort(right);

  output = MergeSort(left, right);

  return true;
}

bool ChetverikovaEShellSortSimpleMergeSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace chetverikova_e_shell_sort_simple_merge
