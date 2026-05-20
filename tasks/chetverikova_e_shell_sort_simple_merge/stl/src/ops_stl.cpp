#include "chetverikova_e_shell_sort_simple_merge/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include "chetverikova_e_shell_sort_simple_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace chetverikova_e_shell_sort_simple_merge {

ChetverikovaEShellSortSimpleMergeSTL::ChetverikovaEShellSortSimpleMergeSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ChetverikovaEShellSortSimpleMergeSTL::ValidationImpl() {
  return !(GetInput().empty());
}

bool ChetverikovaEShellSortSimpleMergeSTL::PreProcessingImpl() {
  return true;
}

void ChetverikovaEShellSortSimpleMergeSTL::ShellSort(std::vector<int> &data) {
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

bool ChetverikovaEShellSortSimpleMergeSTL::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  if (input.empty()) {
    output.clear();
    return true;
  }

  const size_t num_threads = ppc::util::GetNumThreads();
  const size_t n = input.size();

  const size_t block_size = (n + num_threads - 1) / num_threads;
  std::vector<std::vector<int>> blocks(num_threads);
  std::vector<std::thread> threads(num_threads);

  for (size_t tid = 0; tid < num_threads; ++tid) {
    threads[tid] = std::thread([&, tid]() {
      size_t start = tid * block_size;
      size_t end = std::min(start + block_size, n);

      if (start >= n) {
        return;
      }

      std::vector<int> local;
      local.reserve(end - start);

      for (size_t i = start; i < end; ++i) {
        local.push_back(input[i]);
      }

      ShellSort(local);

      blocks[tid] = std::move(local);
    });
  }

  for (auto &th : threads) {
    th.join();
  }

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

bool ChetverikovaEShellSortSimpleMergeSTL::PostProcessingImpl() {
  return true;
}

}  // namespace chetverikova_e_shell_sort_simple_merge
