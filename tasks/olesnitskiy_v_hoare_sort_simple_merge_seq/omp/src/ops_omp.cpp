#include "olesnitskiy_v_hoare_sort_simple_merge_seq/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stack>
#include <utility>
#include <vector>

namespace olesnitskiy_v_hoare_sort_simple_merge_seq {

OlesnitskiyVHoareSortSimpleMergeOMP::OlesnitskiyVHoareSortSimpleMergeOMP(const std::vector<int> &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

int OlesnitskiyVHoareSortSimpleMergeOMP::HoarePartition(std::vector<int> &array, int left, int right) {
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

void OlesnitskiyVHoareSortSimpleMergeOMP::HoareQuickSort(std::vector<int> &array, int left, int right) {
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

void OlesnitskiyVHoareSortSimpleMergeOMP::SimpleMerge(const std::vector<int> &source, std::vector<int> &destination,
                                                      size_t left, size_t middle, size_t right) {
  size_t left_index = left;
  size_t right_index = middle;
  size_t dest_index = left;

  while (left_index < middle && right_index < right) {
    if (source[left_index] <= source[right_index]) {
      destination[dest_index++] = source[left_index++];
    } else {
      destination[dest_index++] = source[right_index++];
    }
  }

  while (left_index < middle) {
    destination[dest_index++] = source[left_index++];
  }

  while (right_index < right) {
    destination[dest_index++] = source[right_index++];
  }
}

void OlesnitskiyVHoareSortSimpleMergeOMP::SortBlocks(std::vector<int> &data, size_t block_size, int num_threads) {
  const size_t size = data.size();
  const size_t block_count = (size + block_size - 1) / block_size;
  const auto block_count_i64 = static_cast<std::int64_t>(block_count);

#pragma omp parallel for default(none) shared(data, size, block_size, block_count_i64) schedule(static) \
    num_threads(num_threads)
  for (std::int64_t block_index = 0; block_index < block_count_i64; ++block_index) {
    size_t block_start = static_cast<size_t>(block_index) * block_size;
    size_t block_end = std::min(block_start + block_size, size);
    if ((block_end - block_start) > 1) {
      HoareQuickSort(data, static_cast<int>(block_start), static_cast<int>(block_end - 1));
    }
  }
}

void OlesnitskiyVHoareSortSimpleMergeOMP::MergeSortedBlocks(std::vector<int> &data, size_t block_size,
                                                            int num_threads) {
  const size_t size = data.size();
  std::vector<int> buffer(size);
  bool data_is_source = true;

  for (size_t merge_width = block_size; merge_width < size; merge_width *= 2) {
    const size_t chunk_width = merge_width * 2;
    const size_t chunk_count = (size + chunk_width - 1) / chunk_width;
    const auto chunk_count_i64 = static_cast<std::int64_t>(chunk_count);
    const std::vector<int> &source = data_is_source ? data : buffer;
    std::vector<int> &destination = data_is_source ? buffer : data;

#pragma omp parallel for default(none) shared(source, destination, size, merge_width, chunk_width, chunk_count_i64) \
    schedule(static) num_threads(num_threads)
    for (std::int64_t chunk_index = 0; chunk_index < chunk_count_i64; ++chunk_index) {
      size_t left = static_cast<size_t>(chunk_index) * chunk_width;
      size_t middle = std::min(left + merge_width, size);
      size_t right = std::min(left + chunk_width, size);

      if (middle < right) {
        SimpleMerge(source, destination, left, middle, right);
        continue;
      }

      std::copy(source.begin() + static_cast<std::ptrdiff_t>(left), source.begin() + static_cast<std::ptrdiff_t>(right),
                destination.begin() + static_cast<std::ptrdiff_t>(left));
    }

    data_is_source = !data_is_source;
    if (merge_width > (size / 2)) {
      break;
    }
  }

  if (!data_is_source) {
    data.swap(buffer);
  }
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

  constexpr size_t kBlockSize = 64;
  const int num_threads = 2;
  SortBlocks(data_, kBlockSize, num_threads);
  MergeSortedBlocks(data_, kBlockSize, num_threads);

  return true;
}

bool OlesnitskiyVHoareSortSimpleMergeOMP::PostProcessingImpl() {
  if (!std::ranges::is_sorted(data_)) {
    return false;
  }
  GetOutput() = data_;
  return true;
}

}  // namespace olesnitskiy_v_hoare_sort_simple_merge_seq
