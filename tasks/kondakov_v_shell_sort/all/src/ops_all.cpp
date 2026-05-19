#include "kondakov_v_shell_sort/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <limits>
#include <thread>
#include <utility>
#include <vector>

#include "kondakov_v_shell_sort/common/include/common.hpp"
#include "util/include/util.hpp"

namespace kondakov_v_shell_sort {

namespace {

size_t CalcPartsCount(size_t data_size, int requested_threads);

template <class F>
void ParallelForIndex(size_t count, int requested_threads, F body) {
  if (count == 0) {
    return;
  }

  const size_t workers_count = std::min(count, static_cast<size_t>(std::max(1, requested_threads)));
  std::atomic_size_t next_idx = 0;
  std::vector<std::thread> workers;
  workers.reserve(workers_count);

  for (size_t worker = 0; worker < workers_count; ++worker) {
    workers.emplace_back([&]() {
      while (true) {
        const size_t idx = next_idx.fetch_add(1);
        if (idx >= count) {
          break;
        }
        body(idx);
      }
    });
  }

  for (auto &worker : workers) {
    worker.join();
  }
}

struct ChunkBounds {
  size_t begin;
  size_t end;
};

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

void MergeManySortedChunks(const std::vector<int> &flat_data, const std::vector<int> &counts,
                           std::vector<int> &result) {
  if (counts.empty()) {
    result.clear();
    return;
  }

  size_t offset = 0;
  std::vector<int> merged;
  for (int count : counts) {
    const size_t chunk_size = static_cast<size_t>(std::max(0, count));
    std::vector<int> chunk(flat_data.begin() + static_cast<std::ptrdiff_t>(offset),
                           flat_data.begin() + static_cast<std::ptrdiff_t>(offset + chunk_size));

    if (merged.empty()) {
      merged = std::move(chunk);
    } else {
      std::vector<int> tmp(merged.size() + chunk.size());
      SimpleMerge(merged, chunk, tmp);
      merged = std::move(tmp);
    }

    offset += chunk_size;
  }

  result = std::move(merged);
}

std::vector<ChunkBounds> BuildBalancedChunks(size_t data_size, size_t chunks_count) {
  std::vector<ChunkBounds> bounds(chunks_count);
  for (size_t i = 0; i < chunks_count; ++i) {
    bounds[i] = ChunkBounds{.begin = (i * data_size) / chunks_count, .end = ((i + 1) * data_size) / chunks_count};
  }
  return bounds;
}

void SortLocalWithStlThreads(std::vector<int> &local_data, int requested_threads) {
  const size_t parts = CalcPartsCount(local_data.size(), requested_threads);
  if (parts <= 1) {
    ShellSort(local_data);
    return;
  }

  const std::vector<ChunkBounds> bounds = BuildBalancedChunks(local_data.size(), parts);
  std::vector<std::vector<int>> runs(parts);
  ParallelForIndex(parts, requested_threads, [&](size_t part) {
    runs[part] = std::vector<int>(local_data.begin() + static_cast<std::ptrdiff_t>(bounds[part].begin),
                                  local_data.begin() + static_cast<std::ptrdiff_t>(bounds[part].end));
    ShellSort(runs[part]);
  });

  std::vector<int> merged = std::move(runs.front());
  for (size_t i = 1; i < runs.size(); ++i) {
    std::vector<int> tmp(merged.size() + runs[i].size());
    SimpleMerge(merged, runs[i], tmp);
    merged = std::move(tmp);
  }
  local_data = std::move(merged);
}

bool IsMpiSizeRepresentable(size_t size) {
  return size <= static_cast<size_t>(std::numeric_limits<int>::max());
}

size_t CalcPartsCount(size_t data_size, int requested_threads) {
  if (data_size == 0) {
    return 0;
  }
  const int threads = std::max(1, requested_threads);
  return std::min(static_cast<size_t>(threads), data_size);
}

}  // namespace

KondakovVShellSortALL::KondakovVShellSortALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool KondakovVShellSortALL::ValidationImpl() {
  return !GetInput().empty();
}

bool KondakovVShellSortALL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool KondakovVShellSortALL::RunImpl() {
  std::vector<int> &data = GetOutput();
  if (data.size() <= 1) {
    return true;
  }

  if (!IsMpiSizeRepresentable(data.size())) {
    return false;
  }

  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int global_size = (rank == 0) ? static_cast<int>(data.size()) : 0;
  MPI_Bcast(&global_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (global_size <= 0) {
    return false;
  }

  std::vector<int> global_data;
  if (rank == 0) {
    global_data = data;
  } else {
    global_data.resize(static_cast<size_t>(global_size));
  }
  MPI_Bcast(global_data.data(), global_size, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> counts(static_cast<size_t>(world_size), 0);
  std::vector<int> displs(static_cast<size_t>(world_size), 0);
  for (int proc = 0; proc < world_size; ++proc) {
    const int begin = (proc * global_size) / world_size;
    const int end = ((proc + 1) * global_size) / world_size;
    counts[static_cast<size_t>(proc)] = end - begin;
    displs[static_cast<size_t>(proc)] = begin;
  }

  std::vector<int> local_data(static_cast<size_t>(counts[static_cast<size_t>(rank)]));
  MPI_Scatterv(global_data.data(), counts.data(), displs.data(), MPI_INT, local_data.data(),
               counts[static_cast<size_t>(rank)], MPI_INT, 0, MPI_COMM_WORLD);

  SortLocalWithStlThreads(local_data, ppc::util::GetNumThreads());

  std::vector<int> gathered;
  if (rank == 0) {
    gathered.resize(static_cast<size_t>(global_size));
  }
  MPI_Gatherv(local_data.data(), counts[static_cast<size_t>(rank)], MPI_INT, gathered.data(), counts.data(),
              displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> sorted_result;
  if (rank == 0) {
    MergeManySortedChunks(gathered, counts, sorted_result);
    if (sorted_result.size() != static_cast<size_t>(global_size)) {
      return false;
    }
  } else {
    sorted_result.resize(static_cast<size_t>(global_size));
  }

  MPI_Bcast(sorted_result.data(), global_size, MPI_INT, 0, MPI_COMM_WORLD);
  data = std::move(sorted_result);

  return std::ranges::is_sorted(data);
}

bool KondakovVShellSortALL::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace kondakov_v_shell_sort
