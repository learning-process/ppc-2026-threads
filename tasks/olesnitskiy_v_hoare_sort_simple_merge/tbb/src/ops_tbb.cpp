#include "olesnitskiy_v_hoare_sort_simple_merge/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>

#include <algorithm>
#include <cstddef>
#include <stack>
#include <utility>
#include <vector>

#include "olesnitskiy_v_hoare_sort_simple_merge/common/include/common.hpp"

namespace olesnitskiy_v_hoare_sort_simple_merge {

namespace {

constexpr size_t kBlockSize = 64;

}  // namespace

OlesnitskiyVHoareSortSimpleMergeTBB::OlesnitskiyVHoareSortSimpleMergeTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

int OlesnitskiyVHoareSortSimpleMergeTBB::HoarePartition(std::vector<int> &array, int left, int right) {
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

void OlesnitskiyVHoareSortSimpleMergeTBB::HoareQuickSort(std::vector<int> &array, int left, int right) {
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

void OlesnitskiyVHoareSortSimpleMergeTBB::SimpleMerge(const std::vector<int> &source, std::vector<int> &destination,
                                                      size_t left, size_t middle, size_t right) {
  size_t left_index = left;
  size_t right_index = middle;
  size_t destination_index = left;

  while (left_index < middle && right_index < right) {
    if (source[left_index] <= source[right_index]) {
      destination[destination_index++] = source[left_index++];
    } else {
      destination[destination_index++] = source[right_index++];
    }
  }

  while (left_index < middle) {
    destination[destination_index++] = source[left_index++];
  }

  while (right_index < right) {
    destination[destination_index++] = source[right_index++];
  }
}

bool OlesnitskiyVHoareSortSimpleMergeTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool OlesnitskiyVHoareSortSimpleMergeTBB::PreProcessingImpl() {
  data_ = GetInput();
  GetOutput().clear();
  return true;
}

bool OlesnitskiyVHoareSortSimpleMergeTBB::RunImpl() {
  if (data_.size() <= 1) {
    return true;
  }

  const size_t size = data_.size();
  const size_t block_count = (size + kBlockSize - 1) / kBlockSize;

  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<size_t>(0, block_count), [this, size](const auto &range) {
    for (size_t block_index = range.begin(); block_index != range.end(); ++block_index) {
      size_t block_start = block_index * kBlockSize;
      size_t block_end = std::min(block_start + kBlockSize, size);
      if ((block_end - block_start) > 1) {
        HoareQuickSort(data_, static_cast<int>(block_start), static_cast<int>(block_end - 1));
      }
    }
  });

  for (size_t merge_width = kBlockSize; merge_width < size; merge_width *= 2) {
    std::vector<int> merged_data(size);
    const size_t merge_count = (size + (2 * merge_width) - 1) / (2 * merge_width);

    oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<size_t>(0, merge_count),
                              [this, size, merge_width, &merged_data](const auto &range) {
      for (size_t merge_index = range.begin(); merge_index != range.end(); ++merge_index) {
        size_t left = merge_index * 2 * merge_width;
        size_t middle = std::min(left + merge_width, size);
        size_t right = std::min(left + (2 * merge_width), size);

        if (middle < right) {
          SimpleMerge(data_, merged_data, left, middle, right);
        } else {
          std::copy(data_.begin() + static_cast<std::ptrdiff_t>(left),
                    data_.begin() + static_cast<std::ptrdiff_t>(right),
                    merged_data.begin() + static_cast<std::ptrdiff_t>(left));
        }
      }
    });

    data_.swap(merged_data);
  }

  return true;
}

bool OlesnitskiyVHoareSortSimpleMergeTBB::PostProcessingImpl() {
  const bool is_sorted = std::ranges::is_sorted(data_);
  if (is_sorted) {
    GetOutput() = data_;
  }
  return is_sorted;
}

}  // namespace olesnitskiy_v_hoare_sort_simple_merge
