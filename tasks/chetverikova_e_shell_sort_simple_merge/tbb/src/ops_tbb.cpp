#include "chetverikova_e_shell_sort_simple_merge/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

#include "chetverikova_e_shell_sort_simple_merge/common/include/common.hpp"

namespace chetverikova_e_shell_sort_simple_merge {

ChetverikovaEShellSortSimpleMergeTBB::ChetverikovaEShellSortSimpleMergeTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ChetverikovaEShellSortSimpleMergeTBB::ValidationImpl() {
  return !(GetInput().empty());
}

bool ChetverikovaEShellSortSimpleMergeTBB::PreProcessingImpl() {
  return true;
}

void ChetverikovaEShellSortSimpleMergeTBB::ShellSort(std::vector<int> &data) {
  if (data.empty()) {
    return;
  }

  size_t n = data.size();
  for (size_t gap = n / 2; gap > 0; gap /= 2) {
    for (size_t i = gap; i < n; i++) {
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

bool ChetverikovaEShellSortSimpleMergeTBB::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  if (input.empty()) {
    output.clear();
    return true;
  }

  const size_t n = input.size();

  const size_t num_threads = std::min<size_t>(4, n);
  const size_t block_size = (n + num_threads - 1) / num_threads;
  std::vector<std::vector<int>> blocks(num_threads);

  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_threads), [&](const tbb::blocked_range<size_t> &r) {
    for (size_t block_id = r.begin(); block_id < r.end(); ++block_id) {
      size_t start = block_id * block_size;
      size_t end = std::min(start + block_size, n);

      if (start >= n) {
        continue;
      }

      std::vector<int> local;
      local.reserve(end - start);

      for (size_t i = start; i < end; ++i) {
        local.push_back(input[i]);
      }

      ShellSort(local);

      blocks[block_id] = std::move(local);
    }
  });

  std::vector<int> result = std::move(blocks[0]);

  for (size_t i = 1; i < num_threads; ++i) {
    if (blocks[i].empty()) {
      continue;
    }

    std::vector<int> tmp(result.size() + blocks[i].size());

    std::merge(result.begin(), result.end(), blocks[i].begin(), blocks[i].end(), tmp.begin());

    result.swap(tmp);
  }

  output = std::move(result);
  return true;
}

bool ChetverikovaEShellSortSimpleMergeTBB::PostProcessingImpl() {
  return true;
}

}  // namespace chetverikova_e_shell_sort_simple_merge
