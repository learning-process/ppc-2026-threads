#include "kondakov_v_shell_sort/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

#include "kondakov_v_shell_sort/common/include/common.hpp"
#include "util/include/util.hpp"

namespace kondakov_v_shell_sort {

namespace {

void ShellSort(std::vector<int> &data) {
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

void SimpleMerge(const std::vector<int> &left, const std::vector<int> &right, std::vector<int> &result) {
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
}

size_t CalcPartsCount(size_t data_size, int requested_threads) {
  if (data_size == 0) {
    return 0;
  }
  const int threads = std::max(1, requested_threads);
  return std::min(static_cast<size_t>(threads), data_size);
}

}  // namespace

KondakovVShellSortOMP::KondakovVShellSortOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool KondakovVShellSortOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool KondakovVShellSortOMP::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool KondakovVShellSortOMP::RunImpl() {
  std::vector<int> &data = GetOutput();
  if (data.size() <= 1) {
    return true;
  }

  const int num_threads = ppc::util::GetNumThreads();
  const size_t parts = CalcPartsCount(data.size(), num_threads);
  if (parts <= 1) {
    ShellSort(data);
    return std::ranges::is_sorted(data);
  }

  std::vector<std::vector<int>> segments(parts);

#pragma omp parallel for default(none) shared(data, segments, parts) num_threads(num_threads) schedule(static)
  for (size_t part = 0; part < parts; ++part) {
    const size_t begin = (part * data.size()) / parts;
    const size_t end = ((part + 1) * data.size()) / parts;
    segments[part] = std::vector<int>(data.begin() + static_cast<std::ptrdiff_t>(begin),
                                      data.begin() + static_cast<std::ptrdiff_t>(end));
    ShellSort(segments[part]);
  }

  std::vector<int> merged = std::move(segments[0]);
  for (size_t i = 1; i < parts; ++i) {
    std::vector<int> tmp(merged.size() + segments[i].size());
    SimpleMerge(merged, segments[i], tmp);
    merged = std::move(tmp);
  }

  data = std::move(merged);
  return std::ranges::is_sorted(data);
}

bool KondakovVShellSortOMP::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace kondakov_v_shell_sort
