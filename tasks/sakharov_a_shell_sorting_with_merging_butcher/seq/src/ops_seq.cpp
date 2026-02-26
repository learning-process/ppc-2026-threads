#include "sakharov_a_shell_sorting_with_merging_butcher/seq/include/ops_seq.hpp"

#include <cstddef>
#include <utility>
#include <vector>

#include "sakharov_a_shell_sorting_with_merging_butcher/common/include/common.hpp"

namespace sakharov_a_shell_sorting_with_merging_butcher {

namespace {

void SplitByParity(const std::vector<int> &source, std::vector<int> &even, std::vector<int> &odd) {
  even.reserve((source.size() + 1) / 2);
  odd.reserve(source.size() / 2);

  for (std::size_t index = 0; index < source.size(); ++index) {
    if (index % 2 == 0) {
      even.push_back(source[index]);
    } else {
      odd.push_back(source[index]);
    }
  }
}

std::vector<int> Interleave(const std::vector<int> &even, const std::vector<int> &odd) {
  std::vector<int> result;
  result.reserve(even.size() + odd.size());

  std::size_t even_index = 0;
  std::size_t odd_index = 0;

  while (even_index < even.size() || odd_index < odd.size()) {
    if (even_index < even.size()) {
      result.push_back(even[even_index]);
      ++even_index;
    }
    if (odd_index < odd.size()) {
      result.push_back(odd[odd_index]);
      ++odd_index;
    }
  }

  return result;
}

void FixOddEvenAdjacents(std::vector<int> &data) {
  for (std::size_t index = 1; index + 1 < data.size(); index += 2) {
    if (data[index] > data[index + 1]) {
      std::swap(data[index], data[index + 1]);
    }
  }
}

}  // namespace

SakharovAShellButcherSEQ::SakharovAShellButcherSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SakharovAShellButcherSEQ::ValidationImpl() {
  return IsValidInput(GetInput());
}

bool SakharovAShellButcherSEQ::PreProcessingImpl() {
  GetOutput().assign(GetInput().size(), 0);
  return true;
}

void SakharovAShellButcherSEQ::ShellSort(std::vector<int> &data) {
  if (data.empty()) {
    return;
  }

  for (std::size_t gap = data.size() / 2; gap > 0; gap /= 2) {
    for (std::size_t i = gap; i < data.size(); ++i) {
      int temp = data[i];
      std::size_t j = i;
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

  std::size_t left_index = 0;
  std::size_t right_index = 0;

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

  SplitByParity(left, left_even, left_odd);
  SplitByParity(right, right_even, right_odd);

  std::vector<int> even_merged = MergeSortedVectors(left_even, right_even);
  std::vector<int> odd_merged = MergeSortedVectors(left_odd, right_odd);

  std::vector<int> result = Interleave(even_merged, odd_merged);
  FixOddEvenAdjacents(result);

  return result;
}

bool SakharovAShellButcherSEQ::RunImpl() {
  const auto &input = GetInput();

  if (input.empty()) {
    GetOutput().clear();
    return true;
  }

  const std::size_t middle = input.size() / 2;

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
