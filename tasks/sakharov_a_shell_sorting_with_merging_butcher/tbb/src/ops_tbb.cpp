#include "sakharov_a_shell_sorting_with_merging_butcher/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <cstddef>
#include <vector>

#include "sakharov_a_shell_sorting_with_merging_butcher/common/include/common.hpp"
#include "util/include/util.hpp"

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
  const auto chunk_count = bounds.size() - 1;

  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, chunk_count), [&](const tbb::blocked_range<std::size_t> &range) {
    for (std::size_t chunk = range.begin(); chunk != range.end(); ++chunk) {
      auto begin = data.begin() + static_cast<std::ptrdiff_t>(bounds[chunk]);
      auto end = data.begin() + static_cast<std::ptrdiff_t>(bounds[chunk + 1]);
      std::sort(begin, end);
    }
  });
}

void ParallelMergePass(const std::vector<int> &source, std::vector<int> &destination,
                       const std::vector<std::size_t> &bounds, std::size_t width) {
  const std::size_t chunk_count = bounds.size() - 1;
  const std::size_t merge_count = (chunk_count + (2 * width) - 1) / (2 * width);

  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, merge_count), [&](const tbb::blocked_range<std::size_t> &range) {
    for (std::size_t merge_index = range.begin(); merge_index != range.end(); ++merge_index) {
      const std::size_t left_chunk = merge_index * 2 * width;
      const std::size_t mid = std::min(left_chunk + width, chunk_count);
      const std::size_t right = std::min(left_chunk + (2 * width), chunk_count);

      const std::size_t begin_index = bounds[left_chunk];
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
  });
}

std::vector<int> ParallelSortAndMerge(const std::vector<int> &input) {
  std::vector<int> source = input;
  const auto bounds = BuildChunkBounds(source.size(), std::max(1, ppc::util::GetNumThreads()));
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

SakharovAShellButcherTBB::SakharovAShellButcherTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool SakharovAShellButcherTBB::ValidationImpl() {
  return IsValidInput(GetInput());
}

bool SakharovAShellButcherTBB::PreProcessingImpl() {
  GetOutput().assign(GetInput().size(), 0);
  return true;
}

bool SakharovAShellButcherTBB::RunImpl() {
  const auto &input = GetInput();

  if (input.empty()) {
    GetOutput().clear();
    return true;
  }

  GetOutput() = ParallelSortAndMerge(input);
  return true;
}

bool SakharovAShellButcherTBB::PostProcessingImpl() {
  return true;
}

}  // namespace sakharov_a_shell_sorting_with_merging_butcher
