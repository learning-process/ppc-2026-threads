#include "ovchinnikov_m_shell_sort_batcher_merge_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

#include "ovchinnikov_m_shell_sort_batcher_merge_seq/common/include/common.hpp"

namespace ovchinnikov_m_shell_sort_batcher_merge_seq {

namespace {

std::size_t NextPowerOfTwo(std::size_t value) {
  std::size_t power = 1;
  while (power < value) {
    power <<= 1;
  }
  return power;
}

std::vector<int> MergeSorted(const std::vector<int> &first, const std::vector<int> &second) {
  std::vector<int> merged;
  merged.reserve(first.size() + second.size());

  std::size_t left_index = 0;
  std::size_t right_index = 0;
  while (left_index < first.size() && right_index < second.size()) {
    if (first[left_index] <= second[right_index]) {
      merged.push_back(first[left_index]);
      ++left_index;
    } else {
      merged.push_back(second[right_index]);
      ++right_index;
    }
  }

  while (left_index < first.size()) {
    merged.push_back(first[left_index]);
    ++left_index;
  }
  while (right_index < second.size()) {
    merged.push_back(second[right_index]);
    ++right_index;
  }

  return merged;
}

bool IsEvenPosition(std::size_t index) {
  return (index % 2) == 0;
}

void SplitByParity(const std::vector<int> &input, std::vector<int> &even, std::vector<int> &odd) {
  even.reserve((input.size() + 1) / 2);
  odd.reserve(input.size() / 2);
  for (std::size_t i = 0; i < input.size(); ++i) {
    if (IsEvenPosition(i)) {
      even.push_back(input[i]);
    } else {
      odd.push_back(input[i]);
    }
  }
}

std::vector<int> InterleaveByParity(const std::vector<int> &even, const std::vector<int> &odd) {
  std::vector<int> merged(even.size() + odd.size());
  std::size_t even_index = 0;
  std::size_t odd_index = 0;
  for (std::size_t i = 0; i < merged.size(); ++i) {
    if (IsEvenPosition(i)) {
      merged[i] = even[even_index];
      ++even_index;
    } else {
      merged[i] = odd[odd_index];
      ++odd_index;
    }
  }
  return merged;
}

std::pair<std::vector<int>, std::vector<int>> SplitInHalf(const std::vector<int> &data) {
  const auto middle = data.begin() + static_cast<std::ptrdiff_t>(data.size() / 2);
  std::vector<int> left(data.begin(), middle);
  std::vector<int> right(middle, data.end());
  return {std::move(left), std::move(right)};
}

void CompareExchangeOddPairs(std::vector<int> &data) {
  for (std::size_t i = 1; i + 1 < data.size(); i += 2) {
    if (data[i] > data[i + 1]) {
      std::swap(data[i], data[i + 1]);
    }
  }
}

}  // namespace

OvchinnikovMShellSortBatcherMergeSEQ::OvchinnikovMShellSortBatcherMergeSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

void OvchinnikovMShellSortBatcherMergeSEQ::ShellSort(std::vector<int> &data) {
  const std::size_t n = data.size();
  for (std::size_t gap = n / 2; gap > 0; gap /= 2) {
    for (std::size_t i = gap; i < n; ++i) {
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

std::vector<int> OvchinnikovMShellSortBatcherMergeSEQ::BatcherOddEvenMerge(const std::vector<int> &left,
                                                                           const std::vector<int> &right) {
  if (left.empty()) {
    return right;
  }
  if (right.empty()) {
    return left;
  }

  if (left.size() != right.size() || left.size() <= 1) {
    return MergeSorted(left, right);
  }

  std::vector<int> left_even;
  std::vector<int> left_odd;
  std::vector<int> right_even;
  std::vector<int> right_odd;
  SplitByParity(left, left_even, left_odd);
  SplitByParity(right, right_even, right_odd);

  std::vector<int> merged_even = MergeSorted(left_even, right_even);
  std::vector<int> merged_odd = MergeSorted(left_odd, right_odd);

  std::vector<int> merged = InterleaveByParity(merged_even, merged_odd);
  CompareExchangeOddPairs(merged);

  return merged;
}

bool OvchinnikovMShellSortBatcherMergeSEQ::ValidationImpl() {
  return true;
}

bool OvchinnikovMShellSortBatcherMergeSEQ::PreProcessingImpl() {
  data_ = GetInput();
  return true;
}

bool OvchinnikovMShellSortBatcherMergeSEQ::RunImpl() {
  constexpr std::size_t kMinSizeToSort = 2;
  if (data_.size() < kMinSizeToSort) {
    return true;
  }

  const std::size_t original_size = data_.size();
  const std::size_t padded_size = NextPowerOfTwo(original_size);
  data_.resize(padded_size, std::numeric_limits<int>::max());

  auto [left, right] = SplitInHalf(data_);

  ShellSort(left);
  ShellSort(right);

  data_ = BatcherOddEvenMerge(left, right);
  data_.resize(original_size);
  return true;
}

bool OvchinnikovMShellSortBatcherMergeSEQ::PostProcessingImpl() {
  if (!std::ranges::is_sorted(data_.begin(), data_.end())) {
    return false;
  }
  GetOutput() = data_;
  return true;
}

}  // namespace ovchinnikov_m_shell_sort_batcher_merge_seq
