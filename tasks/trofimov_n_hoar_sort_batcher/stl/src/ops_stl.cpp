#include "trofimov_n_hoar_sort_batcher/stl/include/ops_stl.hpp"

#include <algorithm>
#include <thread>
#include <vector>

#include "trofimov_n_hoar_sort_batcher/common/include/common.hpp"

namespace trofimov_n_hoar_sort_batcher {

namespace {

unsigned int GetThreadsCount(std::size_t size) {
  unsigned int threads_count = std::thread::hardware_concurrency();
  if (threads_count == 0) {
    threads_count = 2;
  }
  return std::min<unsigned int>(threads_count, static_cast<unsigned int>(size));
}

std::vector<std::pair<std::size_t, std::size_t>> BuildRanges(std::size_t size, unsigned int threads_count) {
  const std::size_t chunk_size = (size + threads_count - 1) / static_cast<std::size_t>(threads_count);
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
  ranges.reserve(threads_count);

  for (std::size_t begin = 0; begin < size; begin += chunk_size) {
    const std::size_t end = std::min(begin + chunk_size, size);
    ranges.emplace_back(begin, end);
  }
  return ranges;
}

void SortRangesInParallel(std::vector<int> &arr, const std::vector<std::pair<std::size_t, std::size_t>> &ranges) {
  std::vector<std::thread> workers;
  workers.reserve(ranges.size());
  for (const auto &[begin, end] : ranges) {
    workers.emplace_back([&arr, begin, end]() { std::sort(arr.begin() + begin, arr.begin() + end); });
  }
  for (auto &worker : workers) {
    worker.join();
  }
}

void MergeSortedRanges(std::vector<int> &arr, std::vector<std::pair<std::size_t, std::size_t>> &ranges) {
  while (ranges.size() > 1) {
    std::vector<std::pair<std::size_t, std::size_t>> merged_ranges;
    merged_ranges.reserve((ranges.size() + 1) / 2);

    for (std::size_t i = 0; i < ranges.size(); i += 2) {
      if (i + 1 >= ranges.size()) {
        merged_ranges.push_back(ranges[i]);
        continue;
      }

      const auto [left_begin, left_end] = ranges[i];
      const auto [right_begin, right_end] = ranges[i + 1];
      std::inplace_merge(arr.begin() + left_begin, arr.begin() + left_end, arr.begin() + right_end);
      merged_ranges.emplace_back(left_begin, right_end);
    }

    ranges.swap(merged_ranges);
  }
}

void ParallelSortByChunks(std::vector<int> &arr) {
  const std::size_t size = arr.size();
  if (size <= 1) {
    return;
  }

  const unsigned int threads_count = GetThreadsCount(size);
  auto ranges = BuildRanges(size, threads_count);
  SortRangesInParallel(arr, ranges);
  MergeSortedRanges(arr, ranges);
}

}  // namespace

TrofimovNHoarSortBatcherSTL::TrofimovNHoarSortBatcherSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TrofimovNHoarSortBatcherSTL::ValidationImpl() {
  return true;
}

bool TrofimovNHoarSortBatcherSTL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool TrofimovNHoarSortBatcherSTL::RunImpl() {
  auto &data = GetOutput();
  ParallelSortByChunks(data);
  return true;
}

bool TrofimovNHoarSortBatcherSTL::PostProcessingImpl() {
  return true;
}

}  // namespace trofimov_n_hoar_sort_batcher
