#include "chetverikova_e_shell_sort_simple_merge/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <utility>
#include <vector>

#include "chetverikova_e_shell_sort_simple_merge/common/include/common.hpp"

namespace chetverikova_e_shell_sort_simple_merge {

ChetverikovaEShellSortSimpleMergeOMP::ChetverikovaEShellSortSimpleMergeOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ChetverikovaEShellSortSimpleMergeOMP::ValidationImpl() {
  return !(GetInput().empty());
}

bool ChetverikovaEShellSortSimpleMergeOMP::PreProcessingImpl() {
  return true;
}

void ChetverikovaEShellSortSimpleMergeOMP::ShellSort(std::vector<int> &data) {
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

bool ChetverikovaEShellSortSimpleMergeOMP::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  if (input.empty()) {
    output.clear();
    return true;
  }

  output = input;
  const std::size_t threads = std::max(1, omp_get_max_threads());
  const std::size_t counts_parts = std::min<std::size_t>(threads, input.size());
  const std::size_t len_part = input.size() / counts_parts;
  const std::size_t rem = input.size() % counts_parts;

  std::vector<size_t> ind_parts;
  ind_parts.push_back(0);

  for (size_t i = 0; i < counts_parts; ++i) {
    ind_parts.push_back(ind_parts.back() + len_part);
    if (i < rem) {
      ind_parts[i + 1]++;
    }
  }
  std::vector<std::vector<int>> local_buffers(counts_parts);

#pragma omp parallel for default(none) shared(input, output, ind_parts, local_buffers, counts_parts) schedule(static)
  for (size_t i = 0; i < counts_parts; ++i) {
    auto left = static_cast<std::ptrdiff_t>(ind_parts[i]);
    auto right = static_cast<std::ptrdiff_t>(ind_parts[i + 1]);

    std::vector<int> temp(input.begin() + left, input.begin() + right);
    ShellSort(temp);
    local_buffers[i] = std::move(temp);
  }

  output = std::move(local_buffers[0]);
  for (size_t i = 1; i < counts_parts; ++i) {
    std::vector<int> merged;
    merged.reserve(output.size() + local_buffers[i].size());
    std::merge(output.begin(), output.end(), local_buffers[i].begin(), local_buffers[i].end(),
               std::back_inserter(merged));
    output = std::move(merged);
  }

  return true;
}

bool ChetverikovaEShellSortSimpleMergeOMP::PostProcessingImpl() {
  return true;
}

}  // namespace chetverikova_e_shell_sort_simple_merge
