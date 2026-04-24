#include "sakharov_a_shell_sorting_with_merging_butcher/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <thread>
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

void ShellSort(std::vector<int>::iterator begin, std::vector<int>::iterator end) {
  const auto size = static_cast<std::size_t>(end - begin);
  if (size <= 1) {
    return;
  }

  for (std::size_t gap = size / 2; gap > 0; gap /= 2) {
    for (std::size_t i = gap; i < size; ++i) {
      const int temp = *(begin + static_cast<std::ptrdiff_t>(i));
      std::size_t j = i;
      while (j >= gap && *(begin + static_cast<std::ptrdiff_t>(j - gap)) > temp) {
        *(begin + static_cast<std::ptrdiff_t>(j)) = *(begin + static_cast<std::ptrdiff_t>(j - gap));
        j -= gap;
      }
      *(begin + static_cast<std::ptrdiff_t>(j)) = temp;
    }
  }
}

void SortChunks(std::vector<int> &data, const std::vector<std::size_t> &bounds) {
  const auto chunk_count = bounds.size() - 1;
  std::vector<std::thread> threads;
  threads.reserve(chunk_count);

  for (std::size_t chunk = 0; chunk < chunk_count; ++chunk) {
    threads.emplace_back([&data, &bounds, chunk]() {
      auto begin = data.begin() + static_cast<std::ptrdiff_t>(bounds[chunk]);
      auto end = data.begin() + static_cast<std::ptrdiff_t>(bounds[chunk + 1]);
      ShellSort(begin, end);
    });
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

void MergePassWorker(const std::vector<int> &source, std::vector<int> &destination,
                     const std::vector<std::size_t> &bounds, std::size_t width, std::size_t begin_merge,
                     std::size_t end_merge) {
  const std::size_t chunk_count = bounds.size() - 1;

  for (std::size_t merge_index = begin_merge; merge_index < end_merge; ++merge_index) {
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
}

void MergePass(const std::vector<int> &source, std::vector<int> &destination, const std::vector<std::size_t> &bounds,
               std::size_t width, int requested_threads) {
  const std::size_t chunk_count = bounds.size() - 1;
  const std::size_t merge_count = (chunk_count + (2 * width) - 1) / (2 * width);
  const auto worker_count =
      std::min<std::size_t>(merge_count, static_cast<std::size_t>(std::max(1, requested_threads)));

  std::vector<std::thread> threads;
  threads.reserve(worker_count);

  const std::size_t base_chunk = merge_count / worker_count;
  const std::size_t remainder = merge_count % worker_count;
  std::size_t current = 0;

  for (std::size_t worker = 0; worker < worker_count; ++worker) {
    const std::size_t count = base_chunk + (worker < remainder ? 1 : 0);
    const std::size_t next = current + count;
    threads.emplace_back(MergePassWorker, std::cref(source), std::ref(destination), std::cref(bounds), width, current,
                         next);
    current = next;
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

std::vector<int> ParallelShellSortAndMerge(const std::vector<int> &input) {
  std::vector<int> source = input;
  const int thread_count = std::max(1, ppc::util::GetNumThreads());
  const auto bounds = BuildChunkBounds(source.size(), thread_count);
  const std::size_t chunk_count = bounds.size() - 1;

  SortChunks(source, bounds);
  if (chunk_count == 1) {
    return source;
  }

  std::vector<int> destination(source.size());
  for (std::size_t width = 1; width < chunk_count; width *= 2) {
    MergePass(source, destination, bounds, width, thread_count);
    source.swap(destination);
  }

  return source;
}

}  // namespace

SakharovAShellButcherSTL::SakharovAShellButcherSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool SakharovAShellButcherSTL::ValidationImpl() {
  return IsValidInput(GetInput());
}

bool SakharovAShellButcherSTL::PreProcessingImpl() {
  GetOutput().assign(GetInput().size(), 0);
  return true;
}

bool SakharovAShellButcherSTL::RunImpl() {
  const auto &input = GetInput();

  if (input.empty()) {
    GetOutput().clear();
    return true;
  }

  GetOutput() = ParallelShellSortAndMerge(input);
  return true;
}

bool SakharovAShellButcherSTL::PostProcessingImpl() {
  return true;
}

}  // namespace sakharov_a_shell_sorting_with_merging_butcher
