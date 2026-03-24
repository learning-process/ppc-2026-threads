#include "olesnitskiy_v_hoare_sort_simple_merge_seq/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stack>
#include <utility>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

#include "util/include/util.hpp"

namespace olesnitskiy_v_hoare_sort_simple_merge_seq {

namespace {

constexpr std::size_t kBlockSize = 64;

int GetActualThreadCount() {
#ifdef _OPENMP
  return std::max(1, std::min(ppc::util::GetNumThreads(), omp_get_max_threads()));
#else
  return 1;
#endif
}

template <typename In, typename Out>
void CopyRange(In first, In last, Out result) {
  std::ranges::copy(first, last, result);
}

}  // namespace

OlesnitskiyVHoareSortSimpleMergeOMP::OlesnitskiyVHoareSortSimpleMergeOMP(const std::vector<int> &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

OlesnitskiyVHoareSortSimpleMergeOMP::~OlesnitskiyVHoareSortSimpleMergeOMP() {
#if defined(_OPENMP) && (_OPENMP >= 201811)
  omp_pause_resource_all(omp_pause_hard);
#endif
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

void OlesnitskiyVHoareSortSimpleMergeOMP::SortBlocks(std::vector<int> &data, size_t block_size, int num_threads) {
  const size_t size = data.size();
  const size_t block_count = (size + block_size - 1) / block_size;
  const auto block_count_i64 = static_cast<std::int64_t>(block_count);

#ifdef _OPENMP
  const int actual_threads = std::max(1, num_threads);
#  pragma omp parallel for default(none) shared(data, size, block_size, block_count_i64) schedule(static) \
      num_threads(actual_threads)
#endif
  for (std::int64_t block_index = 0; block_index < block_count_i64; ++block_index) {
    const size_t block_start = static_cast<size_t>(block_index) * block_size;
    const size_t block_end = std::min(block_start + block_size, size);
    if ((block_end - block_start) > 1) {
      HoareQuickSort(data, static_cast<int>(block_start), static_cast<int>(block_end - 1));
    }
  }
}

void OlesnitskiyVHoareSortSimpleMergeOMP::MergeSortedBlocks(std::vector<int> &data, size_t block_size,
                                                            int num_threads) {
  const size_t size = data.size();
  if (size <= 1) {
    return;
  }

  for (size_t merge_width = block_size; merge_width < size; merge_width *= 2) {
    std::vector<int> merged_data(size);
    const size_t pair_width = merge_width * 2;
    const size_t pair_count = (size + pair_width - 1) / pair_width;
    const auto pair_count_i64 = static_cast<std::int64_t>(pair_count);

#ifdef _OPENMP
    const int actual_threads = std::max(1, num_threads);
#  pragma omp parallel for default(none) shared(data, merged_data, size, merge_width, pair_width, pair_count_i64) \
      schedule(static) num_threads(actual_threads)
#endif
    for (std::int64_t pair_index = 0; pair_index < pair_count_i64; ++pair_index) {
      const size_t left = static_cast<size_t>(pair_index) * pair_width;
      const size_t middle = std::min(left + merge_width, size);
      const size_t right = std::min(left + pair_width, size);

      if (middle < right) {
        const std::vector<int> left_part(data.begin() + static_cast<std::ptrdiff_t>(left),
                                         data.begin() + static_cast<std::ptrdiff_t>(middle));
        const std::vector<int> right_part(data.begin() + static_cast<std::ptrdiff_t>(middle),
                                          data.begin() + static_cast<std::ptrdiff_t>(right));
        const std::vector<int> merged_part = SimpleMerge(left_part, right_part);
        CopyRange(merged_part.begin(), merged_part.end(), merged_data.begin() + static_cast<std::ptrdiff_t>(left));
      } else {
        CopyRange(data.begin() + static_cast<std::ptrdiff_t>(left), data.begin() + static_cast<std::ptrdiff_t>(right),
                  merged_data.begin() + static_cast<std::ptrdiff_t>(left));
      }
    }

    data.swap(merged_data);
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

  const int num_threads = GetActualThreadCount();
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
