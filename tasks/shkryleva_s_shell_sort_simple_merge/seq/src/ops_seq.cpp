#include "shkryleva_s_shell_sort_simple_merge/seq/include/ops_seq.hpp"

#include <cstddef>
#include <vector>

#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

ShkrylevaSShellMergeSEQ::ShkrylevaSShellMergeSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ShkrylevaSShellMergeSEQ::ValidationImpl() {
  return true;
}

bool ShkrylevaSShellMergeSEQ::PreProcessingImpl() {
  input_data_ = GetInput();
  output_data_.clear();
  return true;
}

void ShkrylevaSShellMergeSEQ::ShellSort(std::vector<int> &data) {
  if (data.empty()) {
    return;
  }

  for (size_t gap = data.size() / 2; gap > 0; gap /= 2) {
    for (size_t i = gap; i < data.size(); ++i) {
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

std::vector<int> ShkrylevaSShellMergeSEQ::MergeSortedVectors(const std::vector<int> &left,
                                                             const std::vector<int> &right) {
  std::vector<int> merged;
  merged.reserve(left.size() + right.size());

  size_t left_index = 0;
  size_t right_index = 0;

  while (left_index < left.size() && right_index < right.size()) {
    if (left[left_index] <= right[right_index]) {
      merged.push_back(left[left_index]);
      ++left_index;
    } else {
      merged.push_back(right[right_index]);
      ++right_index;
    }
  }

  while (left_index < left.size()) {
    merged.push_back(left[left_index]);
    ++left_index;
  }

  while (right_index < right.size()) {
    merged.push_back(right[right_index]);
    ++right_index;
  }

  return merged;
}

bool ShkrylevaSShellMergeSEQ::RunImpl() {
  if (input_data_.empty()) {
    output_data_.clear();
    return true;
  }

  const size_t middle = input_data_.size() / 2;

  std::vector<int> left(input_data_.begin(), input_data_.begin() + static_cast<std::ptrdiff_t>(middle));
  std::vector<int> right(input_data_.begin() + static_cast<std::ptrdiff_t>(middle), input_data_.end());

  // 1. Сортируем каждую половину Шеллом
  ShellSort(left);
  ShellSort(right);

  // 2. Обычное линейное слияние
  output_data_ = MergeSortedVectors(left, right);

  return true;
}

bool ShkrylevaSShellMergeSEQ::PostProcessingImpl() {
  GetOutput() = output_data_;
  return true;
}

}  // namespace shkryleva_s_shell_sort_simple_merge
