#include "sakharov_a_shell_sorting_with_merging_butcher/seq/include/ops_seq.hpp"

#include <cstddef>
#include <utility>
#include <vector>

#include "sakharov_a_shell_sorting_with_merging_butcher/common/include/common.hpp"

namespace sakharov_a_shell_sorting_with_merging_butcher {

SakharovAShellButcherSEQ::SakharovAShellButcherSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SakharovAShellButcherSEQ::ValidationImpl() {
  return IsValidInput(GetInput());
}

bool SakharovAShellButcherSEQ::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

void SakharovAShellButcherSEQ::ShellSort(std::vector<int> &data) {
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

std::vector<int> SakharovAShellButcherSEQ::MergeSortedVectors(const std::vector<int> &left,
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

std::vector<int> SakharovAShellButcherSEQ::BatcherOddEvenMerge(const std::vector<int> &left,
                                                               const std::vector<int> &right) {
  std::vector<int> left_even;
  std::vector<int> left_odd;
  std::vector<int> right_even;
  std::vector<int> right_odd;

  left_even.reserve((left.size() + 1) / 2);
  left_odd.reserve(left.size() / 2);
  right_even.reserve((right.size() + 1) / 2);
  right_odd.reserve(right.size() / 2);

  for (size_t index = 0; index < left.size(); ++index) {
    if (index % 2 == 0) {
      left_even.push_back(left[index]);
    } else {
      left_odd.push_back(left[index]);
    }
  }

  for (size_t index = 0; index < right.size(); ++index) {
    if (index % 2 == 0) {
      right_even.push_back(right[index]);
    } else {
      right_odd.push_back(right[index]);
    }
  }

  std::vector<int> even_merged = MergeSortedVectors(left_even, right_even);
  std::vector<int> odd_merged = MergeSortedVectors(left_odd, right_odd);

  std::vector<int> result;
  result.reserve(left.size() + right.size());
  size_t even_index = 0;
  size_t odd_index = 0;

  while (even_index < even_merged.size() || odd_index < odd_merged.size()) {
    if (even_index < even_merged.size()) {
      result.push_back(even_merged[even_index]);
      ++even_index;
    }
    if (odd_index < odd_merged.size()) {
      result.push_back(odd_merged[odd_index]);
      ++odd_index;
    }
  }

  for (size_t index = 1; index + 1 < result.size(); index += 2) {
    if (result[index] > result[index + 1]) {
      std::swap(result[index], result[index + 1]);
    }
  }

  return result;
}

bool SakharovAShellButcherSEQ::RunImpl() {
  const auto &input = GetInput();

  if (input.empty()) {
    GetOutput().clear();
    return true;
  }

  const size_t middle = input.size() / 2;

  std::vector<int> left(input.begin(), input.begin() + static_cast<std::ptrdiff_t>(middle));
  std::vector<int> right(input.begin() + static_cast<std::ptrdiff_t>(middle), input.end());

  ShellSort(left);
  ShellSort(right);

  auto output = BatcherOddEvenMerge(left, right);
  ShellSort(output);
  GetOutput() = std::move(output);
  return true;
}

bool SakharovAShellButcherSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace sakharov_a_shell_sorting_with_merging_butcher
