#include "yushkova_p_hoare_sorting_simple_merging/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stack>
#include <thread>
#include <utility>
#include <vector>

#include "yushkova_p_hoare_sorting_simple_merging/common/include/common.hpp"

namespace yushkova_p_hoare_sorting_simple_merging {

namespace {

constexpr std::size_t kBlockSize = 64;

std::size_t GetThreadCount(std::size_t task_count) {
  if (task_count == 0) {
    return 0;
  }
  const unsigned int hardware_threads = std::thread::hardware_concurrency();
  const std::size_t available_threads = hardware_threads == 0 ? 2 : static_cast<std::size_t>(hardware_threads);
  return std::min(task_count, available_threads);
}

template <class Function>
void RunInThreads(std::size_t task_count, Function function) {
  const std::size_t thread_count = GetThreadCount(task_count);
  if (thread_count <= 1) {
    for (std::size_t task_index = 0; task_index < task_count; ++task_index) {
      function(task_index);
    }
    return;
  }

  std::vector<std::thread> threads;
  threads.reserve(thread_count);
  for (std::size_t thread_index = 0; thread_index < thread_count; ++thread_index) {
    // NOLINTNEXTLINE(bugprone-exception-escape)
    threads.emplace_back([thread_index, thread_count, task_count, &function]() noexcept {
      for (std::size_t task_index = thread_index; task_index < task_count; task_index += thread_count) {
        function(task_index);
      }
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }
}

std::vector<int> MakeIntVector(const std::vector<std::size_t> &values) {
  std::vector<int> result(values.size());
  for (std::size_t i = 0; i < values.size(); ++i) {
    result[i] = static_cast<int>(values[i]);
  }
  return result;
}

void BuildDistribution(std::size_t total_size, int mpi_size, std::vector<std::size_t> &chunk_sizes,
                       std::vector<std::size_t> &offsets) {
  const std::size_t base_size = total_size / static_cast<std::size_t>(mpi_size);
  const std::size_t remainder = total_size % static_cast<std::size_t>(mpi_size);
  for (int rank = 0; rank < mpi_size; ++rank) {
    const auto index = static_cast<std::size_t>(rank);
    chunk_sizes[index] = base_size + (index < remainder ? 1U : 0U);
    offsets[index] = rank == 0 ? 0 : offsets[index - 1] + chunk_sizes[index - 1];
  }
}

}  // namespace

YushkovaPHoareSortingSimpleMergingALL::YushkovaPHoareSortingSimpleMergingALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

int YushkovaPHoareSortingSimpleMergingALL::HoarePartition(std::vector<int> &values, int left, int right) {
  const int pivot = values[left + ((right - left) / 2)];
  int i = left - 1;
  int j = right + 1;

  while (true) {
    ++i;
    while (values[i] < pivot) {
      ++i;
    }

    --j;
    while (values[j] > pivot) {
      --j;
    }

    if (i >= j) {
      return j;
    }

    std::swap(values[i], values[j]);
  }
}

void YushkovaPHoareSortingSimpleMergingALL::HoareQuickSort(std::vector<int> &values, int left, int right) {
  std::stack<std::pair<int, int>> ranges;
  ranges.emplace(left, right);

  while (!ranges.empty()) {
    auto [current_left, current_right] = ranges.top();
    ranges.pop();

    if (current_left >= current_right) {
      continue;
    }

    const int partition_index = HoarePartition(values, current_left, current_right);

    if ((partition_index - current_left) > (current_right - (partition_index + 1))) {
      ranges.emplace(current_left, partition_index);
      ranges.emplace(partition_index + 1, current_right);
    } else {
      ranges.emplace(partition_index + 1, current_right);
      ranges.emplace(current_left, partition_index);
    }
  }
}

void YushkovaPHoareSortingSimpleMergingALL::SimpleMerge(const std::vector<int> &source, std::vector<int> &destination,
                                                        std::size_t left, std::size_t middle, std::size_t right) {
  std::size_t left_index = left;
  std::size_t right_index = middle;
  std::size_t destination_index = left;

  while (left_index < middle && right_index < right) {
    if (source[left_index] <= source[right_index]) {
      destination[destination_index++] = source[left_index++];
    } else {
      destination[destination_index++] = source[right_index++];
    }
  }

  while (left_index < middle) {
    destination[destination_index++] = source[left_index++];
  }

  while (right_index < right) {
    destination[destination_index++] = source[right_index++];
  }
}

void YushkovaPHoareSortingSimpleMergingALL::SortLocalStlParallel(std::vector<int> &values) {
  if (values.size() <= 1) {
    return;
  }

  const std::size_t size = values.size();
  const std::size_t block_count = (size + kBlockSize - 1) / kBlockSize;

  RunInThreads(block_count, [&values, size](std::size_t block_index) {
    const std::size_t block_start = block_index * kBlockSize;
    const std::size_t block_end = std::min(block_start + kBlockSize, size);
    if ((block_end - block_start) > 1) {
      HoareQuickSort(values, static_cast<int>(block_start), static_cast<int>(block_end - 1));
    }
  });

  for (std::size_t merge_width = kBlockSize; merge_width < size; merge_width *= 2) {
    std::vector<int> merged_data(size);
    const std::size_t merge_count = (size + (2 * merge_width) - 1) / (2 * merge_width);

    RunInThreads(merge_count, [&values, size, merge_width, &merged_data](std::size_t merge_index) {
      const std::size_t left = merge_index * 2 * merge_width;
      const std::size_t middle = std::min(left + merge_width, size);
      const std::size_t right = std::min(left + (2 * merge_width), size);

      if (middle < right) {
        SimpleMerge(values, merged_data, left, middle, right);
      } else {
        std::copy(values.begin() + static_cast<std::ptrdiff_t>(left),
                  values.begin() + static_cast<std::ptrdiff_t>(right),
                  merged_data.begin() + static_cast<std::ptrdiff_t>(left));
      }
    });

    values.swap(merged_data);
  }
}

void YushkovaPHoareSortingSimpleMergingALL::MergeGatheredChunks(std::vector<int> &values,
                                                                const std::vector<std::size_t> &chunk_sizes,
                                                                const std::vector<std::size_t> &offsets) {
  std::vector<int> merged_data(values.size());
  for (std::size_t rank = 1; rank < chunk_sizes.size(); ++rank) {
    const std::size_t left = 0;
    const std::size_t middle = offsets[rank];
    const std::size_t right = middle + chunk_sizes[rank];
    SimpleMerge(values, merged_data, left, middle, right);
    std::copy(merged_data.begin(), merged_data.begin() + static_cast<std::ptrdiff_t>(right), values.begin());
  }
}

void YushkovaPHoareSortingSimpleMergingALL::BroadcastVector(std::vector<int> &values, int rank) {
  auto size = static_cast<std::uint64_t>(values.size());
  MPI_Bcast(&size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  if (rank != 0) {
    values.resize(static_cast<std::size_t>(size));
  }
  if (size > 0) {
    MPI_Bcast(values.data(), static_cast<int>(size), MPI_INT, 0, MPI_COMM_WORLD);
  }
}

bool YushkovaPHoareSortingSimpleMergingALL::ValidationImpl() {
  return !GetInput().empty();
}

bool YushkovaPHoareSortingSimpleMergingALL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool YushkovaPHoareSortingSimpleMergingALL::RunImpl() {
  int rank = 0;
  int mpi_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  const std::size_t total_size = GetOutput().size();
  if (total_size <= 1) {
    return true;
  }

  std::vector<std::size_t> chunk_sizes(static_cast<std::size_t>(mpi_size));
  std::vector<std::size_t> offsets(static_cast<std::size_t>(mpi_size));
  BuildDistribution(total_size, mpi_size, chunk_sizes, offsets);

  const std::vector<int> send_counts = MakeIntVector(chunk_sizes);
  const std::vector<int> send_offsets = MakeIntVector(offsets);

  std::vector<int> local_data(chunk_sizes[static_cast<std::size_t>(rank)]);
  MPI_Scatterv(rank == 0 ? GetOutput().data() : nullptr, send_counts.data(), send_offsets.data(), MPI_INT,
               local_data.data(), send_counts[static_cast<std::size_t>(rank)], MPI_INT, 0, MPI_COMM_WORLD);

  SortLocalStlParallel(local_data);

  std::vector<int> gathered_data;
  if (rank == 0) {
    gathered_data.resize(total_size);
  }
  MPI_Gatherv(local_data.data(), static_cast<int>(local_data.size()), MPI_INT,
              rank == 0 ? gathered_data.data() : nullptr, send_counts.data(), send_offsets.data(), MPI_INT, 0,
              MPI_COMM_WORLD);

  if (rank == 0) {
    MergeGatheredChunks(gathered_data, chunk_sizes, offsets);
    GetOutput() = std::move(gathered_data);
  }

  BroadcastVector(GetOutput(), rank);
  return std::ranges::is_sorted(GetOutput());
}

bool YushkovaPHoareSortingSimpleMergingALL::PostProcessingImpl() {
  return !GetOutput().empty() && std::ranges::is_sorted(GetOutput());
}

}  // namespace yushkova_p_hoare_sorting_simple_merging
