#include "ovchinnikov_m_shell_sort_batcher_merge/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

#include "ovchinnikov_m_shell_sort_batcher_merge/common/include/common.hpp"

namespace ovchinnikov_m_shell_sort_batcher_merge {

namespace {

std::size_t NextPowerOfTwo(std::size_t value) {
  std::size_t power = 1;
  while (power < value) {
    power <<= 1;
  }
  return power;
}

std::vector<int> MergeSorted(const std::vector<int> &first, const std::vector<int> &second) {
  std::vector<int> merged;
  merged.reserve(first.size() + second.size());

  std::size_t left_index = 0;
  std::size_t right_index = 0;
  while (left_index < first.size() && right_index < second.size()) {
    if (first[left_index] <= second[right_index]) {
      merged.push_back(first[left_index]);
      ++left_index;
    } else {
      merged.push_back(second[right_index]);
      ++right_index;
    }
  }

  while (left_index < first.size()) {
    merged.push_back(first[left_index]);
    ++left_index;
  }
  while (right_index < second.size()) {
    merged.push_back(second[right_index]);
    ++right_index;
  }

  return merged;
}

bool IsEvenPosition(std::size_t index) {
  return (index % 2) == 0;
}

void SplitByParity(const std::vector<int> &input, std::vector<int> &even, std::vector<int> &odd) {
  even.reserve((input.size() + 1) / 2);
  odd.reserve(input.size() / 2);
  for (std::size_t i = 0; i < input.size(); ++i) {
    if (IsEvenPosition(i)) {
      even.push_back(input[i]);
    } else {
      odd.push_back(input[i]);
    }
  }
}

std::vector<int> InterleaveByParity(const std::vector<int> &even, const std::vector<int> &odd) {
  std::vector<int> merged(even.size() + odd.size());
  std::size_t even_index = 0;
  std::size_t odd_index = 0;
  for (std::size_t i = 0; i < merged.size(); ++i) {
    if (IsEvenPosition(i)) {
      merged[i] = even[even_index];
      ++even_index;
    } else {
      merged[i] = odd[odd_index];
      ++odd_index;
    }
  }
  return merged;
}

void CompareExchangeOddPairs(std::vector<int> &data) {
  for (std::size_t i = 1; i + 1 < data.size(); i += 2) {
    if (data[i] > data[i + 1]) {
      std::swap(data[i], data[i + 1]);
    }
  }
}

int LargestPowerOfTwoNotGreaterThan(int value) {
  int power = 1;
  while ((power * 2) <= value) {
    power *= 2;
  }
  return power;
}

}  // namespace

OvchinnikovMShellSortBatcherMergeALL::OvchinnikovMShellSortBatcherMergeALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

void OvchinnikovMShellSortBatcherMergeALL::ShellSort(std::vector<int> &data) {
  const std::size_t n = data.size();
  for (std::size_t gap = n / 2; gap > 0; gap /= 2) {
    for (std::size_t i = gap; i < n; ++i) {
      const int temp = data[i];
      std::size_t j = i;
      while (j >= gap && data[j - gap] > temp) {
        data[j] = data[j - gap];
        j -= gap;
      }
      data[j] = temp;
    }
  }
}

std::vector<int> OvchinnikovMShellSortBatcherMergeALL::BatcherOddEvenMerge(const std::vector<int> &left,
                                                                           const std::vector<int> &right) {
  if (left.empty()) {
    return right;
  }
  if (right.empty()) {
    return left;
  }

  if (left.size() != right.size() || left.size() <= 1) {
    return MergeSorted(left, right);
  }

  std::vector<int> left_even;
  std::vector<int> left_odd;
  std::vector<int> right_even;
  std::vector<int> right_odd;
  SplitByParity(left, left_even, left_odd);
  SplitByParity(right, right_even, right_odd);

  std::vector<int> merged_even = MergeSorted(left_even, right_even);
  std::vector<int> merged_odd = MergeSorted(left_odd, right_odd);

  std::vector<int> merged = InterleaveByParity(merged_even, merged_odd);
  CompareExchangeOddPairs(merged);
  return merged;
}

void OvchinnikovMShellSortBatcherMergeALL::LocalSort(std::vector<int> &local_data) {
  constexpr std::size_t kMinSizeToSort = 2;
  if (local_data.size() < kMinSizeToSort) {
    return;
  }

  const std::size_t half_size = local_data.size() / 2;
  const auto middle = local_data.begin() + static_cast<std::ptrdiff_t>(half_size);
  std::vector<int> left(local_data.begin(), middle);
  std::vector<int> right(middle, local_data.end());

#pragma omp parallel sections default(none) shared(left, right)
  {
#pragma omp section
    {
      ShellSort(left);
    }
#pragma omp section
    {
      ShellSort(right);
    }
  }

  local_data = BatcherOddEvenMerge(left, right);
}

void OvchinnikovMShellSortBatcherMergeALL::TreeMerge(int rank, int active_processes, int chunk_size,
                                                     std::vector<int> &local_data) {
  for (int step = 1; step < active_processes; step *= 2) {
    if ((rank % (2 * step)) == 0) {
      const int sender = rank + step;
      if (sender >= active_processes) {
        continue;
      }

      std::vector<int> received_data(chunk_size * step);
      MPI_Recv(received_data.data(), static_cast<int>(received_data.size()), MPI_INT, sender, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      local_data = BatcherOddEvenMerge(local_data, received_data);
    } else {
      const int receiver = rank - step;
      MPI_Send(local_data.data(), static_cast<int>(local_data.size()), MPI_INT, receiver, 0, MPI_COMM_WORLD);
      break;
    }
  }
}

void OvchinnikovMShellSortBatcherMergeALL::ScatterData(int rank, int active_processes, int chunk_size,
                                                       std::size_t padded_size, std::vector<int> &local_data) {
  if (rank == 0) {
    std::vector<int> padded_data = GetInput();
    padded_data.resize(padded_size, std::numeric_limits<int>::max());

    const auto chunk_size_diff = static_cast<std::ptrdiff_t>(chunk_size);
    std::copy(padded_data.begin(), padded_data.begin() + chunk_size_diff, local_data.begin());

    for (int process = 1; process < active_processes; ++process) {
      const auto offset = static_cast<std::ptrdiff_t>(process * chunk_size);
      MPI_Send(padded_data.data() + offset, chunk_size, MPI_INT, process, 0, MPI_COMM_WORLD);
    }
  } else {
    MPI_Recv(local_data.data(), chunk_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}

bool OvchinnikovMShellSortBatcherMergeALL::ValidationImpl() {
  int is_valid = 1;
  MPI_Bcast(&is_valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return is_valid == 1;
}

bool OvchinnikovMShellSortBatcherMergeALL::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool OvchinnikovMShellSortBatcherMergeALL::RunImpl() {
  int rank = 0;
  int process_count = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &process_count);

  int original_size_int = 0;
  if (rank == 0) {
    original_size_int = static_cast<int>(GetInput().size());
  }
  MPI_Bcast(&original_size_int, 1, MPI_INT, 0, MPI_COMM_WORLD);

  const auto original_size = static_cast<std::size_t>(original_size_int);
  if (original_size <= 1) {
    std::vector<int> result(original_size);
    if (rank == 0 && original_size == 1) {
      result[0] = GetInput()[0];
    }
    if (original_size == 1) {
      MPI_Bcast(result.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    GetOutput() = result;
    return true;
  }

  const int active_processes = LargestPowerOfTwoNotGreaterThan(process_count);
  std::size_t padded_size = NextPowerOfTwo(original_size);
  while (padded_size < static_cast<std::size_t>(active_processes)) {
    padded_size <<= 1;
  }

  const int chunk_size = static_cast<int>(padded_size / static_cast<std::size_t>(active_processes));
  std::vector<int> local_data;

  if (rank < active_processes) {
    local_data.resize(chunk_size);
    ScatterData(rank, active_processes, chunk_size, padded_size, local_data);
    LocalSort(local_data);
    TreeMerge(rank, active_processes, chunk_size, local_data);
  }

  std::vector<int> result(original_size);
  if (rank == 0) {
    std::copy(local_data.begin(), local_data.begin() + static_cast<std::ptrdiff_t>(original_size), result.begin());
  }
  MPI_Bcast(result.data(), original_size_int, MPI_INT, 0, MPI_COMM_WORLD);
  GetOutput() = result;

  return true;
}

bool OvchinnikovMShellSortBatcherMergeALL::PostProcessingImpl() {
  return GetOutput().size() == GetInput().size() && std::ranges::is_sorted(GetOutput().begin(), GetOutput().end());
}

}  // namespace ovchinnikov_m_shell_sort_batcher_merge
