#include "shelenkova_m_shell_sort_simple_merge/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

#include "shelenkova_m_shell_sort_simple_merge/common/include/common.hpp"

namespace shelenkova_m_shell_sort_simple_merge {

namespace {

void ShellSort(std::vector<int> &data) {
  const size_t n = data.size();
  if (n <= 1) {
    return;
  }

  for (size_t gap = n / 2; gap > 0; gap /= 2) {
    for (size_t i = gap; i < n; ++i) {
      int temp = data[i];
      size_t j = i;
      while (j >= gap && data[j - gap] > temp) {
        data[j] = data[j - gap];
        j -= gap;
      }
      data[j] = temp;
    }
  }
}

}  // namespace

ShelenkovaMShellSortSimpleMergeTBB::ShelenkovaMShellSortSimpleMergeTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool ShelenkovaMShellSortSimpleMergeTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool ShelenkovaMShellSortSimpleMergeTBB::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool ShelenkovaMShellSortSimpleMergeTBB::RunImpl() {
  const std::vector<int> &input = GetInput();
  std::vector<int> &output = GetOutput();

  const size_t n = input.size();
  if (n <= 1) {
    return true;
  }

  const int raw_concurrency = tbb::this_task_arena::max_concurrency();
  const size_t num_threads = std::min(n, static_cast<size_t>(std::max(1, raw_concurrency)));
  const size_t block_size = (n + num_threads - 1) / num_threads;

  std::vector<std::vector<int>> blocks(num_threads);

  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_threads), [&](const tbb::blocked_range<size_t> &r) {
    for (size_t block_id = r.begin(); block_id < r.end(); ++block_id) {
      const size_t start = block_id * block_size;
      if (start >= n) {
        continue;
      }
      const size_t end = std::min(start + block_size, n);
      blocks[block_id].assign(input.begin() + static_cast<std::ptrdiff_t>(start),
                              input.begin() + static_cast<std::ptrdiff_t>(end));
      ShellSort(blocks[block_id]);
    }
  });

  std::vector<int> result;
  result.reserve(n);
  for (size_t i = 0; i < num_threads; ++i) {
    if (blocks[i].empty()) {
      continue;
    }
    std::vector<int> tmp(result.size() + blocks[i].size());
    std::merge(result.begin(), result.end(), blocks[i].begin(), blocks[i].end(), tmp.begin());
    result = std::move(tmp);
  }

  output = std::move(result);
  return std::ranges::is_sorted(output);
}

bool ShelenkovaMShellSortSimpleMergeTBB::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace shelenkova_m_shell_sort_simple_merge
