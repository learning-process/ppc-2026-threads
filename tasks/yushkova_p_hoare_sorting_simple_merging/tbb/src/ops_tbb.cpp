#include "yushkova_p_hoare_sorting_simple_merging/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/global_control.h>
#include <oneapi/tbb/parallel_for.h>

#include <algorithm>
#include <stack>
#include <utility>
#include <vector>

#include "util/include/util.hpp"
#include "yushkova_p_hoare_sorting_simple_merging/common/include/common.hpp"

namespace yushkova_p_hoare_sorting_simple_merging {

YushkovaPHoareSortingSimpleMergingTBB::YushkovaPHoareSortingSimpleMergingTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
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
  std::stack<std::pair<int, int>> ranges;
  ranges.emplace(left, right);

  while (!ranges.empty()) {
    auto [current_left, current_right] = ranges.top();
    ranges.pop();

    if (current_left >= current_right) {
      continue;
    }

    const int partition_index = HoarePartition(values, current_left, current_right);

    if ((partition_index - current_left) > (current_right - (partition_index + 1))) {
      ranges.emplace(current_left, partition_index);
      ranges.emplace(partition_index + 1, current_right);
    } else {
      ranges.emplace(partition_index + 1, current_right);
      ranges.emplace(current_left, partition_index);
    }
  }
}

void YushkovaPHoareSortingSimpleMergingTBB::Merge(std::vector<int> &values, int left, int mid, int right) {
  std::vector<int> merged;
  const int merged_size = (right - left) + 1;
  merged.reserve(static_cast<std::size_t>(merged_size));

  int left_index = left;
  int right_index = mid + 1;

  while (left_index <= mid && right_index <= right) {
    if (values[left_index] <= values[right_index]) {
      merged.push_back(values[left_index++]);
    } else {
      merged.push_back(values[right_index++]);
    }
  }

  while (left_index <= mid) {
    merged.push_back(values[left_index++]);
  }

  while (right_index <= right) {
    merged.push_back(values[right_index++]);
  }

  for (std::size_t idx = 0; idx < merged.size(); ++idx) {
    values[static_cast<std::size_t>(left) + idx] = merged[idx];
  }
}

bool YushkovaPHoareSortingSimpleMergingTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool YushkovaPHoareSortingSimpleMergingTBB::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool YushkovaPHoareSortingSimpleMergingTBB::RunImpl() {
  std::vector<int> &values = GetOutput();
  const int n = static_cast<int>(values.size());
  if (n <= 1) {
    return true;
  }

  int num_threads = ppc::util::GetNumThreads();
  if (num_threads < 1) {
    num_threads = 1;
  }

  const int chunks = std::min(num_threads, n);
  if (chunks == 1) {
    HoareQuickSort(values, 0, n - 1);
    return std::ranges::is_sorted(values);
  }

  oneapi::tbb::global_control control(oneapi::tbb::global_control::max_allowed_parallelism,
                                      static_cast<std::size_t>(num_threads));

  std::vector<int> borders(static_cast<std::size_t>(chunks + 1));
  for (int i = 0; i <= chunks; ++i) {
    borders[static_cast<std::size_t>(i)] = (i * n) / chunks;
  }

  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<int>(0, chunks), [&](const oneapi::tbb::blocked_range<int> &r) {
    for (int chunk = r.begin(); chunk < r.end(); ++chunk) {
      const int left = borders[static_cast<std::size_t>(chunk)];
      const int right = borders[static_cast<std::size_t>(chunk) + 1] - 1;
      if (left < right) {
        HoareQuickSort(values, left, right);
      }
    }
  });

  for (int i = 0; i < chunks - 1; ++i) {
    const int mid = borders[static_cast<std::size_t>(i) + 1] - 1;
    const int right = borders[static_cast<std::size_t>(i) + 2] - 1;
    Merge(values, 0, mid, right);
  }

  return std::ranges::is_sorted(values);
}

bool YushkovaPHoareSortingSimpleMergingTBB::PostProcessingImpl() {
  return !GetOutput().empty() && std::ranges::is_sorted(GetOutput());
}

}  // namespace yushkova_p_hoare_sorting_simple_merging
