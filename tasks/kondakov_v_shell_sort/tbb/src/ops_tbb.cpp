#include "kondakov_v_shell_sort/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

#include "kondakov_v_shell_sort/common/include/common.hpp"
#include "util/include/util.hpp"

namespace kondakov_v_shell_sort {

namespace {

struct RunBounds {
  size_t begin;
  size_t end;
};

void SortRunWithShellGaps(std::vector<int> &data) {
  const size_t n = data.size();
  for (size_t gap = n / 2; gap > 0; gap /= 2) {
    for (size_t i = gap; i < n; ++i) {
      int value = data[i];
      size_t j = i;
      while (j >= gap && data[j - gap] > value) {
        data[j] = data[j - gap];
        j -= gap;
      }
      data[j] = value;
    }
  }
}

std::vector<int> MergeTwoSortedRuns(const std::vector<int> &left, const std::vector<int> &right) {
  std::vector<int> result(left.size() + right.size());
  size_t i = 0;
  size_t j = 0;
  size_t k = 0;

  while (i < left.size() && j < right.size()) {
    if (left[i] <= right[j]) {
      result[k++] = left[i++];
    } else {
      result[k++] = right[j++];
    }
  }

  while (i < left.size()) {
    result[k++] = left[i++];
  }

  while (j < right.size()) {
    result[k++] = right[j++];
  }

  return result;
}

size_t GetRunsCount(size_t data_size, int requested_threads) {
  if (data_size == 0) {
    return 0;
  }
  const int threads = std::max(1, requested_threads);
  return std::min(static_cast<size_t>(threads), data_size);
}

std::vector<RunBounds> BuildBalancedRuns(size_t data_size, size_t runs_count) {
  std::vector<RunBounds> bounds(runs_count);
  for (size_t run = 0; run < runs_count; ++run) {
    bounds[run] = RunBounds{.begin = (run * data_size) / runs_count, .end = ((run + 1) * data_size) / runs_count};
  }
  return bounds;
}

std::vector<std::vector<int>> MakeSortedRuns(const std::vector<int> &data, const std::vector<RunBounds> &bounds) {
  std::vector<std::vector<int>> runs(bounds.size());

  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<size_t>(0, bounds.size()),
                            [&](const oneapi::tbb::blocked_range<size_t> &range) {
    for (size_t run = range.begin(); run != range.end(); ++run) {
      const RunBounds &current = bounds[run];
      runs[run] = std::vector<int>(data.begin() + static_cast<std::ptrdiff_t>(current.begin),
                                   data.begin() + static_cast<std::ptrdiff_t>(current.end));
      SortRunWithShellGaps(runs[run]);
    }
  });

  return runs;
}

std::vector<int> MergeRunsByLevels(std::vector<std::vector<int>> runs) {
  while (runs.size() > 1) {
    const size_t pairs_count = runs.size() / 2;
    std::vector<std::vector<int>> next_level(pairs_count + (runs.size() % 2));

    oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<size_t>(0, pairs_count),
                              [&](const oneapi::tbb::blocked_range<size_t> &range) {
      for (size_t pair = range.begin(); pair != range.end(); ++pair) {
        next_level[pair] = MergeTwoSortedRuns(runs[2 * pair], runs[(2 * pair) + 1]);
      }
    });

    if ((runs.size() % 2) != 0) {
      next_level.back() = std::move(runs.back());
    }
    runs = std::move(next_level);
  }

  return std::move(runs.front());
}

}  // namespace

KondakovVShellSortTBB::KondakovVShellSortTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool KondakovVShellSortTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool KondakovVShellSortTBB::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool KondakovVShellSortTBB::RunImpl() {
  std::vector<int> &data = GetOutput();
  if (data.size() <= 1) {
    return true;
  }

  const int num_threads = ppc::util::GetNumThreads();
  const size_t runs_count = GetRunsCount(data.size(), num_threads);
  if (runs_count <= 1) {
    SortRunWithShellGaps(data);
    return std::ranges::is_sorted(data);
  }

  const std::vector<RunBounds> bounds = BuildBalancedRuns(data.size(), runs_count);
  data = MergeRunsByLevels(MakeSortedRuns(data, bounds));
  return std::ranges::is_sorted(data);
}

bool KondakovVShellSortTBB::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace kondakov_v_shell_sort
