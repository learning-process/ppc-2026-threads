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
  GetOutput() = {};
}

int OlesnitskiyVHoareSortSimpleMergeSEQ::HoarePartition(std::vector<int> &array, int left, int right) {
  int pivot = array[left + ((right - left) / 2)];
  int i = left - 1;
  int j = right + 1;

  while (true) {
    i++;
    while (array[i] < pivot) {
      i++;
    }

    j--;
    while (array[j] > pivot) {
      j--;
    }

    if (i >= j) {
      return j;
    }

    std::swap(array[i], array[j]);
  }
}

void OlesnitskiyVHoareSortSimpleMergeSEQ::HoareQuickSort(std::vector<int> &array, int left, int right) {
  std::stack<std::pair<int, int>> stack;
  stack.emplace(left, right);

  while (!stack.empty()) {
    auto [current_left, current_right] = stack.top();
    stack.pop();

    if (current_left >= current_right) {
      continue;
    }

    int middle = HoarePartition(array, current_left, current_right);

    if ((middle - current_left) > (current_right - (middle + 1))) {
      stack.emplace(current_left, middle);
      stack.emplace(middle + 1, current_right);
    } else {
      stack.emplace(middle + 1, current_right);
      stack.emplace(current_left, middle);
    }
  }
}

std::vector<int> OlesnitskiyVHoareSortSimpleMergeSEQ::SimpleMerge(const std::vector<int> &left_part,
                                                                   const std::vector<int> &right_part) {
  std::vector<int> result;
  result.reserve(left_part.size() + right_part.size());

  size_t left_index = 0;
  size_t right_index = 0;

  while (left_index < left_part.size() && right_index < right_part.size()) {
    if (left_part[left_index] <= right_part[right_index]) {
      result.push_back(left_part[left_index]);
      left_index++;
    } else {
      result.push_back(right_part[right_index]);
      right_index++;
    }
  }

  while (left_index < left_part.size()) {
    result.push_back(left_part[left_index]);
    left_index++;
  }

  while (right_index < right_part.size()) {
    result.push_back(right_part[right_index]);
    right_index++;
  }

  return result;
}

bool OlesnitskiyVHoareSortSimpleMergeSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool OlesnitskiyVHoareSortSimpleMergeSEQ::PreProcessingImpl() {
  data_ = GetInput();
  GetOutput().clear();
  return true;
}

bool OlesnitskiyVHoareSortSimpleMergeSEQ::RunImpl() {
  if (data_.size() <= 1) {
    return true;
  }

  constexpr size_t kBlockSize = 64;
  const size_t size = data_.size();

  for (size_t block_start = 0; block_start < size; block_start += kBlockSize) {
    size_t block_end = std::min(block_start + kBlockSize, size);
    if ((block_end - block_start) > 1) {
      HoareQuickSort(data_, static_cast<int>(block_start), static_cast<int>(block_end - 1));
    }
  }

  for (size_t merge_width = kBlockSize; merge_width < size; merge_width *= 2) {
    std::vector<int> merged_data(size);

    for (size_t left = 0; left < size; left += (2 * merge_width)) {
      size_t middle = std::min(left + merge_width, size);
      size_t right = std::min(left + (2 * merge_width), size);

      if (middle < right) {
        std::vector<int> left_part(data_.begin() + static_cast<std::ptrdiff_t>(left),
                                   data_.begin() + static_cast<std::ptrdiff_t>(middle));
        std::vector<int> right_part(data_.begin() + static_cast<std::ptrdiff_t>(middle),
                                    data_.begin() + static_cast<std::ptrdiff_t>(right));
        std::vector<int> merged_part = SimpleMerge(left_part, right_part);
        std::ranges::copy(merged_part, merged_data.begin() + static_cast<std::ptrdiff_t>(left));
      } else {
        std::ranges::copy(data_.begin() + static_cast<std::ptrdiff_t>(left), data_.begin() + static_cast<std::ptrdiff_t>(right),
                          merged_data.begin() + static_cast<std::ptrdiff_t>(left));
      }
    }

    data_.swap(merged_data);
  }

  return true;
}

bool OlesnitskiyVHoareSortSimpleMergeSEQ::PostProcessingImpl() {
  if (!std::ranges::is_sorted(data_)) {
    return false;
  }
  GetOutput() = data_;
  return true;
}

}  //  namespace olesnitskiy_v_hoare_sort_simple_merge_seq
