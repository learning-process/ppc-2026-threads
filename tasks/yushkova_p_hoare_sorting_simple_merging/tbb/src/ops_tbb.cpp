#include "yushkova_p_hoare_sorting_simple_merging/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

#include "yushkova_p_hoare_sorting_simple_merging/common/include/common.hpp"

namespace yushkova_p_hoare_sorting_simple_merging {

namespace {
constexpr size_t kBlockSize = 64;
}  // namespace

YushkovaPHoareSortingSimpleMergingTBB::YushkovaPHoareSortingSimpleMergingTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

int YushkovaPHoareSortingSimpleMergingTBB::HoarePartition(std::vector<int> &values, int left, int right) {
  const int pivot = values[left + ((right - left) / 2)];
  int i = left - 1;
  int j = right + 1;

  while (true) {
    ++i;
    while (values[i] < pivot) {
      ++i;
    }

    --j;
    while (values[j] > pivot) {
      --j;
    }

    if (i >= j) {
      return j;
    }

    std::swap(values[i], values[j]);
  }
}

void YushkovaPHoareSortingSimpleMergingTBB::HoareQuickSort(std::vector<int> &values, int left, int right) {
  if (left >= right || static_cast<size_t>(left) >= values.size() || static_cast<size_t>(right) >= values.size()) {
    return;
  }

  std::vector<std::pair<int, int>> stack;
  stack.emplace_back(left, right);

  while (!stack.empty()) {
    auto [current_left, current_right] = stack.back();
    stack.pop_back();

    if (current_left >= current_right) {
      continue;
    }

    if (current_left < 0 || static_cast<size_t>(current_right) >= values.size()) {
      continue;
    }

    const int partition_index = HoarePartition(values, current_left, current_right);
    if (partition_index < current_left || partition_index > current_right) {
      continue;
    }

    if ((partition_index - current_left) > (current_right - (partition_index + 1))) {
      stack.emplace_back(current_left, partition_index);
      stack.emplace_back(partition_index + 1, current_right);
    } else {
      stack.emplace_back(partition_index + 1, current_right);
      stack.emplace_back(current_left, partition_index);
    }
  }
}

void YushkovaPHoareSortingSimpleMergingTBB::SimpleMerge(const std::vector<int> &source, std::vector<int> &destination,
                                                        size_t left, size_t middle, size_t right) {
  if (left >= right || middle > right || left >= middle) {
    return;
  }

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

void YushkovaPHoareSortingSimpleMergingTBB::SortBlockIfNeeded(std::vector<int> &data, size_t size, size_t block_index) {
  const size_t block_start = block_index * kBlockSize;
  const size_t block_end = std::min(block_start + kBlockSize, size);
  if (block_end > block_start + 1U && block_end <= size && block_start < size) {
    HoareQuickSort(data, static_cast<int>(block_start), static_cast<int>(block_end - 1U));
  }
}

void YushkovaPHoareSortingSimpleMergingTBB::SortBlocks(std::vector<int> &data, size_t size) {
  const size_t block_count = (size + kBlockSize - 1U) / kBlockSize;
  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<size_t>(0U, block_count),
                            [&data, size](const oneapi::tbb::blocked_range<size_t> &range) {
    for (size_t block_index = range.begin(); block_index != range.end(); ++block_index) {
      SortBlockIfNeeded(data, size, block_index);
    }
  });
}

void YushkovaPHoareSortingSimpleMergingTBB::MergeChunk(const std::vector<int> &source, std::vector<int> &destination,
                                                       size_t size, size_t merge_width, size_t merge_index) {
  const size_t left = merge_index * 2U * merge_width;
  const size_t middle = std::min(left + merge_width, size);
  const size_t right = std::min(left + (2U * merge_width), size);

  if (left >= size) {
    return;
  }

  if (middle < right) {
    SimpleMerge(source, destination, left, middle, right);
  } else if (left < right) {
    for (size_t i = left; i < right; ++i) {
      destination[i] = source[i];
    }
  }
}

void YushkovaPHoareSortingSimpleMergingTBB::MergePass(std::vector<int> &data, size_t size, size_t merge_width) {
  std::vector<int> merged_data(size);
  const size_t merge_count = (size + (2U * merge_width) - 1U) / (2U * merge_width);

  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<size_t>(0U, merge_count),
                            [&data, size, merge_width, &merged_data](const oneapi::tbb::blocked_range<size_t> &range) {
    for (size_t merge_index = range.begin(); merge_index != range.end(); ++merge_index) {
      MergeChunk(data, merged_data, size, merge_width, merge_index);
    }
  });

  data.swap(merged_data);
}

bool YushkovaPHoareSortingSimpleMergingTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool YushkovaPHoareSortingSimpleMergingTBB::PreProcessingImpl() {
  data_ = GetInput();
  GetOutput().clear();
  return true;
}

bool YushkovaPHoareSortingSimpleMergingTBB::RunImpl() {
  const size_t size = data_.size();
<<<<<<< HEAD

  if (size == 0) {
    GetOutput().clear();
    return true;
  }

  if (size == 1) {
    GetOutput().assign(1, data_[0]);
    return true;
  }

  SortBlocks(data_, size);
  for (size_t merge_width = kBlockSize; merge_width < size; merge_width *= 2) {
    MergePass(data_, size, merge_width);
  }

  if (std::ranges::is_sorted(data_)) {
=======
  if (size <= 1) {
>>>>>>> 68f721f8 (TBB)
    GetOutput() = data_;
    return true;
  }

<<<<<<< HEAD
=======
  const size_t block_count = (size + kBlockSize - 1) / kBlockSize;

  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<size_t>(0, block_count),
                            [this, size](const oneapi::tbb::blocked_range<size_t> &range) {
                              for (size_t block_index = range.begin(); block_index != range.end(); ++block_index) {
                                const size_t block_start = block_index * kBlockSize;
                                const size_t block_end = std::min(block_start + kBlockSize, size);
                                if (block_end > block_start + 1) {
                                  HoareQuickSort(data_, static_cast<int>(block_start),
                                                 static_cast<int>(block_end - 1));
                                }
                              }
                            });

  std::vector<int> merged_data(size);

  for (size_t merge_width = kBlockSize; merge_width < size; merge_width *= 2) {
    const size_t merge_count = (size + (2 * merge_width) - 1) / (2 * merge_width);

    oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<size_t>(0, merge_count),
                              [this, size, merge_width, &merged_data](const oneapi::tbb::blocked_range<size_t> &range) {
                                for (size_t merge_index = range.begin(); merge_index != range.end(); ++merge_index) {
                                  const size_t left = merge_index * 2 * merge_width;
                                  const size_t middle = std::min(left + merge_width, size);
                                  const size_t right = std::min(left + 2 * merge_width, size);
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

  if (std::is_sorted(data_.begin(), data_.end())) {
    GetOutput().assign(data_.begin(), data_.end());
    return true;
  }
  
>>>>>>> 68f721f8 (TBB)
  return false;
}

bool YushkovaPHoareSortingSimpleMergingTBB::PostProcessingImpl() {
  if (GetOutput().empty()) {
    return false;
  }
  return std::ranges::is_sorted(GetOutput());
}

}  // namespace yushkova_p_hoare_sorting_simple_merging
