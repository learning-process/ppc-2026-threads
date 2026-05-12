#include "olesnitskiy_v_hoare_sort_simple_merge/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <stack>
#include <thread>
#include <utility>
#include <vector>

#include "olesnitskiy_v_hoare_sort_simple_merge/common/include/common.hpp"

namespace olesnitskiy_v_hoare_sort_simple_merge {

namespace {

constexpr size_t kBlockSize = 64;

size_t GetThreadCount(size_t task_count) {
  if (task_count == 0) {
    return 0;
  }

  const unsigned int hardware_threads = std::thread::hardware_concurrency();
  const size_t available_threads = hardware_threads == 0 ? 2 : static_cast<size_t>(hardware_threads);
  return std::min(task_count, available_threads);
}

template <class Function>
void RunInThreads(size_t task_count, Function function) {
  const size_t thread_count = GetThreadCount(task_count);
  if (thread_count <= 1) {
    for (size_t task_index = 0; task_index < task_count; ++task_index) {
      function(task_index);
    }
    return;
  }

  std::vector<std::thread> threads;
  threads.reserve(thread_count);

  for (size_t thread_index = 0; thread_index < thread_count; ++thread_index) {
    threads.emplace_back([thread_index, thread_count, task_count, &function]() {
      for (size_t task_index = thread_index; task_index < task_count; task_index += thread_count) {
        function(task_index);
      }
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }
}

}  // namespace

OlesnitskiyVHoareSortSimpleMergeSTL::OlesnitskiyVHoareSortSimpleMergeSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

int OlesnitskiyVHoareSortSimpleMergeSTL::HoarePartition(std::vector<int> &array, int left, int right) {
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

void OlesnitskiyVHoareSortSimpleMergeSTL::HoareQuickSort(std::vector<int> &array, int left, int right) {
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

void OlesnitskiyVHoareSortSimpleMergeSTL::SimpleMerge(const std::vector<int> &source, std::vector<int> &destination,
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

bool OlesnitskiyVHoareSortSimpleMergeSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool OlesnitskiyVHoareSortSimpleMergeSTL::PreProcessingImpl() {
  data_ = GetInput();
  GetOutput().clear();
  return true;
}

bool OlesnitskiyVHoareSortSimpleMergeSTL::RunImpl() {
  if (data_.size() <= 1) {
    return true;
  }

  const size_t size = data_.size();
  const size_t block_count = (size + kBlockSize - 1) / kBlockSize;

  RunInThreads(block_count, [this, size](size_t block_index) {
    size_t block_start = block_index * kBlockSize;
    size_t block_end = std::min(block_start + kBlockSize, size);
    if ((block_end - block_start) > 1) {
      HoareQuickSort(data_, static_cast<int>(block_start), static_cast<int>(block_end - 1));
    }
  });

  for (size_t merge_width = kBlockSize; merge_width < size; merge_width *= 2) {
    std::vector<int> merged_data(size);
    const size_t merge_count = (size + (2 * merge_width) - 1) / (2 * merge_width);

    RunInThreads(merge_count, [this, size, merge_width, &merged_data](size_t merge_index) {
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
    });

    data_.swap(merged_data);
  }

  return true;
}

bool OlesnitskiyVHoareSortSimpleMergeSTL::PostProcessingImpl() {
  const bool is_sorted = std::ranges::is_sorted(data_);
  if (is_sorted) {
    GetOutput() = data_;
  }
  return is_sorted;
}

}  // namespace olesnitskiy_v_hoare_sort_simple_merge
