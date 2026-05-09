#include "sakharov_a_shell_sorting_with_merging_butcher/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <thread>
#include <utility>
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

std::vector<std::size_t> BuildBoundsFromCounts(const std::vector<int> &counts) {
  std::vector<std::size_t> bounds;
  bounds.reserve(counts.size() + 1);
  bounds.push_back(0);
  for (int count : counts) {
    bounds.push_back(bounds.back() + static_cast<std::size_t>(count));
  }
  return bounds;
}

std::vector<int> BuildMpiChunkCounts(int data_size, int process_count) {
  std::vector<int> counts(static_cast<std::size_t>(process_count), 0);
  const int base = data_size / process_count;
  const int remainder = data_size % process_count;

  for (int rank = 0; rank < process_count; ++rank) {
    counts[static_cast<std::size_t>(rank)] = base + (rank < remainder ? 1 : 0);
  }

  return counts;
}

std::vector<int> BuildDisplacements(const std::vector<int> &counts) {
  std::vector<int> displacements(counts.size(), 0);
  for (std::size_t i = 1; i < counts.size(); ++i) {
    displacements[i] = displacements[i - 1] + counts[i - 1];
  }
  return displacements;
}

void SortChunksOpenMP(std::vector<int> &data, const std::vector<std::size_t> &bounds, int thread_count) {
  const int chunk_count = static_cast<int>(bounds.size()) - 1;

#pragma omp parallel for default(none) shared(data, bounds, chunk_count) num_threads(thread_count) schedule(static)
  for (int chunk = 0; chunk < chunk_count; ++chunk) {
    const auto index = static_cast<std::size_t>(chunk);
    auto begin = data.begin() + static_cast<std::ptrdiff_t>(bounds[index]);
    auto end = data.begin() + static_cast<std::ptrdiff_t>(bounds[index + 1]);
    std::sort(begin, end);
  }
}

void MergeRange(const std::vector<int> &source, std::vector<int> &destination, const std::vector<std::size_t> &bounds,
                std::size_t width, std::size_t merge_index) {
  const std::size_t chunk_count = bounds.size() - 1;
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

void MergePassTBB(const std::vector<int> &source, std::vector<int> &destination, const std::vector<std::size_t> &bounds,
                  std::size_t width) {
  const std::size_t chunk_count = bounds.size() - 1;
  const std::size_t merge_count = (chunk_count + (2 * width) - 1) / (2 * width);

  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, merge_count), [&](const tbb::blocked_range<std::size_t> &range) {
    for (std::size_t merge_index = range.begin(); merge_index != range.end(); ++merge_index) {
      MergeRange(source, destination, bounds, width, merge_index);
    }
  });
}

void MergePassThreadWorker(const std::vector<int> &source, std::vector<int> &destination,
                           const std::vector<std::size_t> &bounds, std::size_t width, std::size_t begin_merge,
                           std::size_t end_merge) {
  for (std::size_t merge_index = begin_merge; merge_index < end_merge; ++merge_index) {
    MergeRange(source, destination, bounds, width, merge_index);
  }
}

void MergePassSTLThreads(const std::vector<int> &source, std::vector<int> &destination,
                         const std::vector<std::size_t> &bounds, std::size_t width, int requested_threads) {
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
    threads.emplace_back(MergePassThreadWorker, std::cref(source), std::ref(destination), std::cref(bounds), width,
                         current, next);
    current = next;
  }

  for (auto &thread : threads) {
    thread.join();
  }
}

std::vector<int> MergeChunksTBB(std::vector<int> source, const std::vector<std::size_t> &bounds) {
  const std::size_t chunk_count = bounds.size() - 1;
  if (chunk_count <= 1) {
    return source;
  }

  std::vector<int> destination(source.size());
  for (std::size_t width = 1; width < chunk_count; width *= 2) {
    MergePassTBB(source, destination, bounds, width);
    source.swap(destination);
  }

  return source;
}

std::vector<int> MergeChunksSTLThreads(std::vector<int> source, const std::vector<std::size_t> &bounds,
                                       int thread_count) {
  const std::size_t chunk_count = bounds.size() - 1;
  if (chunk_count <= 1) {
    return source;
  }

  std::vector<int> destination(source.size());
  for (std::size_t width = 1; width < chunk_count; width *= 2) {
    MergePassSTLThreads(source, destination, bounds, width, thread_count);
    source.swap(destination);
  }

  return source;
}

std::vector<int> SortLocalPart(std::vector<int> data) {
  if (data.empty()) {
    return data;
  }

  const int thread_count = std::max(1, ppc::util::GetNumThreads());
  const auto bounds = BuildChunkBounds(data.size(), thread_count);

  SortChunksOpenMP(data, bounds, thread_count);
  return MergeChunksTBB(std::move(data), bounds);
}

struct MpiRootData {
  int input_size{0};
  std::vector<int> counts;
  std::vector<int> displacements;
};

const int *RootBuffer(const std::vector<int> &buffer, int rank) {
  if (rank != 0) {
    return nullptr;
  }
  return buffer.data();
}

const int *RootInputBuffer(const std::vector<int> &input, int rank) {
  if (rank != 0 || input.empty()) {
    return nullptr;
  }
  return input.data();
}

int *BufferOrNull(std::vector<int> &buffer) {
  if (buffer.empty()) {
    return nullptr;
  }
  return buffer.data();
}

MpiRootData BuildRootData(const std::vector<int> &input, int rank, int process_count) {
  MpiRootData root_data;
  if (rank == 0) {
    root_data.input_size = static_cast<int>(input.size());
    root_data.counts = BuildMpiChunkCounts(root_data.input_size, process_count);
    root_data.displacements = BuildDisplacements(root_data.counts);
  }
  return root_data;
}

std::vector<int> ScatterInput(const std::vector<int> &input, const MpiRootData &root_data, int rank) {
  int local_size = 0;
  MPI_Scatter(RootBuffer(root_data.counts, rank), 1, MPI_INT, &local_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> local_data(static_cast<std::size_t>(local_size));
  MPI_Scatterv(RootInputBuffer(input, rank), RootBuffer(root_data.counts, rank),
               RootBuffer(root_data.displacements, rank), MPI_INT, BufferOrNull(local_data), local_size, MPI_INT, 0,
               MPI_COMM_WORLD);

  return local_data;
}

std::vector<int> GatherSortedData(std::vector<int> &local_data, const MpiRootData &root_data, int rank) {
  std::vector<int> gathered_data;
  if (rank == 0) {
    gathered_data.resize(static_cast<std::size_t>(root_data.input_size));
  }

  const auto local_size = static_cast<int>(local_data.size());
  MPI_Gatherv(BufferOrNull(local_data), local_size, MPI_INT, BufferOrNull(gathered_data),
              RootBuffer(root_data.counts, rank), RootBuffer(root_data.displacements, rank), MPI_INT, 0,
              MPI_COMM_WORLD);

  return gathered_data;
}

std::vector<int> MergeRootData(std::vector<int> gathered_data, const MpiRootData &root_data, int rank) {
  if (rank != 0) {
    return {};
  }

  const auto bounds = BuildBoundsFromCounts(root_data.counts);
  return MergeChunksSTLThreads(std::move(gathered_data), bounds, std::max(1, ppc::util::GetNumThreads()));
}

std::vector<int> BroadcastResult(std::vector<int> result, int rank) {
  int result_size = 0;
  if (rank == 0) {
    result_size = static_cast<int>(result.size());
  }
  MPI_Bcast(&result_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  result.resize(static_cast<std::size_t>(result_size));
  MPI_Bcast(BufferOrNull(result), result_size, MPI_INT, 0, MPI_COMM_WORLD);

  return result;
}

std::vector<int> RunMpiSort(const std::vector<int> &input, int rank, int process_count) {
  const auto root_data = BuildRootData(input, rank, process_count);
  auto local_data = ScatterInput(input, root_data, rank);
  local_data = SortLocalPart(std::move(local_data));
  auto gathered_data = GatherSortedData(local_data, root_data, rank);
  auto result = MergeRootData(std::move(gathered_data), root_data, rank);
  return BroadcastResult(std::move(result), rank);
}

}  // namespace

SakharovAShellButcherALL::SakharovAShellButcherALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool SakharovAShellButcherALL::ValidationImpl() {
  return IsValidInput(GetInput());
}

bool SakharovAShellButcherALL::PreProcessingImpl() {
  GetOutput().assign(GetInput().size(), 0);
  return true;
}

bool SakharovAShellButcherALL::RunImpl() {
  int rank = 0;
  int process_count = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &process_count);

  GetOutput() = RunMpiSort(GetInput(), rank, process_count);
  return true;
}

bool SakharovAShellButcherALL::PostProcessingImpl() {
  return true;
}

}  // namespace sakharov_a_shell_sorting_with_merging_butcher
