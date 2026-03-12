#include "sakharov_a_shell_sorting_with_merging_butcher/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <vector>

#include "sakharov_a_shell_sorting_with_merging_butcher/common/include/common.hpp"

namespace sakharov_a_shell_sorting_with_merging_butcher {

namespace {

constexpr std::size_t kMinParallelChunkSize = 1U << 14;

std::vector<std::size_t> BuildChunkBounds(std::size_t size, int requested_chunks) {
  if (size == 0) {
    return {0};
  }

  const std::size_t max_chunks_by_size = std::max<std::size_t>(1, size / kMinParallelChunkSize);
  const int chunks = std::max(1, std::min<int>(requested_chunks, static_cast<int>(max_chunks_by_size)));

  std::vector<std::size_t> bounds;
  bounds.reserve(static_cast<std::size_t>(chunks) + 1);
  const std::size_t base_chunk = size / static_cast<std::size_t>(chunks);
  const std::size_t remainder = size % static_cast<std::size_t>(chunks);
  const auto chunk_count = static_cast<std::size_t>(chunks);

  bounds.push_back(0);
  for (std::size_t chunk = 0; chunk < chunk_count; ++chunk) {
    const std::size_t chunk_size = base_chunk + (chunk < remainder ? 1 : 0);
    bounds.push_back(bounds.back() + chunk_size);
  }

  return bounds;
}

void ParallelSortChunks(std::vector<int> &data, const std::vector<std::size_t> &bounds) {
  const int chunk_count = static_cast<int>(bounds.size()) - 1;

#pragma omp parallel for default(none) shared(data, bounds, chunk_count) schedule(static)
  for (int chunk = 0; chunk < chunk_count; ++chunk) {
    auto begin = data.begin() + static_cast<std::ptrdiff_t>(bounds[static_cast<std::size_t>(chunk)]);
    auto end = data.begin() + static_cast<std::ptrdiff_t>(bounds[static_cast<std::size_t>(chunk) + 1]);
    std::sort(begin, end);
  }
}

void ParallelMergePass(const std::vector<int> &source, std::vector<int> &destination,
                       const std::vector<std::size_t> &bounds, std::size_t width) {
  const std::size_t chunk_count = bounds.size() - 1;

#pragma omp parallel for default(none) shared(source, destination, bounds, width, chunk_count) schedule(static)
  for (std::size_t left_chunk = 0; left_chunk < chunk_count; left_chunk += 2 * width) {
    auto left = left_chunk;
    const std::size_t mid = std::min(left + width, chunk_count);
    const std::size_t right = std::min(left + (2 * width), chunk_count);

    const std::size_t begin_index = bounds[left];
    const std::size_t middle_index = bounds[mid];
    const std::size_t end_index = bounds[right];

    auto output_begin = destination.begin() + static_cast<std::ptrdiff_t>(begin_index);
    if (mid == right) {
      std::copy(source.begin() + static_cast<std::ptrdiff_t>(begin_index),
                source.begin() + static_cast<std::ptrdiff_t>(end_index), output_begin);
    } else {
      std::merge(source.begin() + static_cast<std::ptrdiff_t>(begin_index),
                 source.begin() + static_cast<std::ptrdiff_t>(middle_index),
                 source.begin() + static_cast<std::ptrdiff_t>(middle_index),
                 source.begin() + static_cast<std::ptrdiff_t>(end_index), output_begin);
    }
  }
}

std::vector<int> ParallelSortAndMerge(const std::vector<int> &input) {
  std::vector<int> source = input;
  const auto bounds = BuildChunkBounds(source.size(), omp_get_max_threads());
  const std::size_t chunk_count = bounds.size() - 1;

  ParallelSortChunks(source, bounds);
  if (chunk_count == 1) {
    return source;
  }

  std::vector<int> destination(source.size());
  for (std::size_t width = 1; width < chunk_count; width *= 2) {
    ParallelMergePass(source, destination, bounds, width);
    source.swap(destination);
  }

  return source;
}

}  // namespace

SakharovAShellButcherOMP::SakharovAShellButcherOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool SakharovAShellButcherOMP::ValidationImpl() {
  return IsValidInput(GetInput());
}

bool SakharovAShellButcherOMP::PreProcessingImpl() {
  GetOutput().assign(GetInput().size(), 0);
  return true;
}

bool SakharovAShellButcherOMP::RunImpl() {
  const auto &input = GetInput();

  if (input.empty()) {
    GetOutput().clear();
    return true;
  }

  GetOutput() = ParallelSortAndMerge(input);
  return true;
}

bool SakharovAShellButcherOMP::PostProcessingImpl() {
  return true;
}

}  // namespace sakharov_a_shell_sorting_with_merging_butcher
