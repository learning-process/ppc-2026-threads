#include "olesnitskiy_v_hoare_sort_simple_merge/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <stack>
#include <utility>
#include <vector>

#include "olesnitskiy_v_hoare_sort_simple_merge/common/include/common.hpp"

namespace olesnitskiy_v_hoare_sort_simple_merge {

OlesnitskiyVHoareSortSimpleMergeOMP::OlesnitskiyVHoareSortSimpleMergeOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

int OlesnitskiyVHoareSortSimpleMergeOMP::HoarePartition(std::vector<int> &array, int left, int right) {
  const int pivot = array[left + ((right - left) / 2)];
  int i = left - 1;
  int j = right + 1;

  while (true) {
    ++i;
    while (array[i] < pivot) {
      ++i;
    }

    --j;
    while (array[j] > pivot) {
      --j;
    }

    if (i >= j) {
      return j;
    }

    std::swap(array[i], array[j]);
  }
}

void OlesnitskiyVHoareSortSimpleMergeOMP::HoareQuickSort(std::vector<int> &array, int left, int right) {
  std::stack<std::pair<int, int>> stack;
  stack.emplace(left, right);

  while (!stack.empty()) {
    auto [current_left, current_right] = stack.top();
    stack.pop();

    if (current_left >= current_right) {
      continue;
    }

    const int middle = HoarePartition(array, current_left, current_right);

    if ((middle - current_left) > (current_right - (middle + 1))) {
      stack.emplace(current_left, middle);
      stack.emplace(middle + 1, current_right);
    } else {
      stack.emplace(middle + 1, current_right);
      stack.emplace(current_left, middle);
    }
  }
}

std::vector<int> OlesnitskiyVHoareSortSimpleMergeOMP::SimpleMerge(const std::vector<int> &left_part,
                                                                  const std::vector<int> &right_part) {
  std::vector<int> result;
  result.reserve(left_part.size() + right_part.size());

  std::size_t left_index = 0;
  std::size_t right_index = 0;

  while (left_index < left_part.size() && right_index < right_part.size()) {
    if (left_part[left_index] <= right_part[right_index]) {
      result.push_back(left_part[left_index]);
      ++left_index;
    } else {
      result.push_back(right_part[right_index]);
      ++right_index;
    }
  }

  while (left_index < left_part.size()) {
    result.push_back(left_part[left_index]);
    ++left_index;
  }

  while (right_index < right_part.size()) {
    result.push_back(right_part[right_index]);
    ++right_index;
  }

  return result;
}

bool OlesnitskiyVHoareSortSimpleMergeOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool OlesnitskiyVHoareSortSimpleMergeOMP::PreProcessingImpl() {
  data_ = GetInput();
  GetOutput().clear();
  return true;
}

bool OlesnitskiyVHoareSortSimpleMergeOMP::RunImpl() {
  if (data_.size() <= 1) {
    return true;
  }

  constexpr std::size_t kBlockSize = 64;
  auto &data = data_;
  const std::size_t size = data.size();

#pragma omp parallel for default(none) shared(size, data)
  for (std::size_t block_start = 0; block_start < size; block_start += kBlockSize) {
    const std::size_t block_end = std::min(block_start + kBlockSize, size);
    if ((block_end - block_start) > 1) {
      HoareQuickSort(data, static_cast<int>(block_start), static_cast<int>(block_end - 1));
    }
  }

  for (std::size_t merge_width = kBlockSize; merge_width < size; merge_width *= 2) {
    std::vector<int> merged_data(size);

#pragma omp parallel for default(none) shared(merge_width, size, merged_data, data)
    for (std::size_t left = 0; left < size; left += (2 * merge_width)) {
      const std::size_t middle = std::min(left + merge_width, size);
      const std::size_t right = std::min(left + (2 * merge_width), size);

      if (middle < right) {
        std::vector<int> left_part(data.begin() + static_cast<std::ptrdiff_t>(left),
                                   data.begin() + static_cast<std::ptrdiff_t>(middle));
        std::vector<int> right_part(data.begin() + static_cast<std::ptrdiff_t>(middle),
                                    data.begin() + static_cast<std::ptrdiff_t>(right));
        std::vector<int> merged_part = SimpleMerge(left_part, right_part);
        for (std::size_t idx = 0; idx < merged_part.size(); ++idx) {
          merged_data[left + idx] = merged_part[idx];
        }
      } else {
        std::copy(data.begin() + static_cast<std::ptrdiff_t>(left), data.begin() + static_cast<std::ptrdiff_t>(right),
                  merged_data.begin() + static_cast<std::ptrdiff_t>(left));
      }
    }

    data.swap(merged_data);
  }

  return true;
}

bool OlesnitskiyVHoareSortSimpleMergeOMP::PostProcessingImpl() {
  if (!std::ranges::is_sorted(data_)) {
    return false;
  }
  GetOutput() = data_;
  return true;
}

}  // namespace olesnitskiy_v_hoare_sort_simple_merge
