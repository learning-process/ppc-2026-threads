#include "tochilin_e_hoar_sort_sim_mer/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <utility>
#include <vector>

#include "tochilin_e_hoar_sort_sim_mer/common/include/common.hpp"
#include "util/include/util.hpp"

namespace tochilin_e_hoar_sort_sim_mer {

namespace {

std::vector<int> BuildCounts(int total_size, int proc_count) {
  std::vector<int> counts(static_cast<std::size_t>(proc_count), 0);
  for (int rank = 0; rank < proc_count; ++rank) {
    counts[static_cast<std::size_t>(rank)] =
        (((rank + 1) * total_size) / proc_count) - ((rank * total_size) / proc_count);
  }
  return counts;
}

std::vector<int> BuildDisplacements(const std::vector<int> &counts) {
  std::vector<int> displs(counts.size(), 0);
  for (std::size_t i = 1; i < counts.size(); ++i) {
    displs[i] = displs[i - 1] + counts[i - 1];
  }
  return displs;
}

std::vector<int> MergeLocalVectors(const std::vector<int> &a, const std::vector<int> &b) {
  std::vector<int> result;
  result.reserve(a.size() + b.size());
  std::ranges::merge(a, b, std::back_inserter(result));
  return result;
}

std::vector<int> MergeAcrossRanks(std::vector<int> local_data, int rank, int proc_count) {
  for (int step = 1; step < proc_count; step *= 2) {
    if ((rank % (step * 2)) == step) {
      const int target_rank = rank - step;
      const int local_size = static_cast<int>(local_data.size());
      MPI_Send(&local_size, 1, MPI_INT, target_rank, 0, MPI_COMM_WORLD);
      if (local_size > 0) {
        MPI_Send(local_data.data(), local_size, MPI_INT, target_rank, 1, MPI_COMM_WORLD);
      }
      break;
    }

    if ((rank % (step * 2)) != 0) {
      continue;
    }

    const int source_rank = rank + step;
    if (source_rank >= proc_count) {
      continue;
    }

    int remote_size = 0;
    MPI_Recv(&remote_size, 1, MPI_INT, source_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    std::vector<int> remote_data(static_cast<std::size_t>(remote_size), 0);
    if (remote_size > 0) {
      MPI_Recv(remote_data.data(), remote_size, MPI_INT, source_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    local_data = MergeLocalVectors(local_data, remote_data);
  }

  return local_data;
}

}  // namespace

TochilinEHoarSortSimMerALL::TochilinEHoarSortSimMerALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TochilinEHoarSortSimMerALL::ValidationImpl() {
  return !GetInput().empty();
}

bool TochilinEHoarSortSimMerALL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

std::pair<int, int> TochilinEHoarSortSimMerALL::Partition(std::vector<int> &arr, int l, int r) {
  int i = l;
  int j = r;
  const int pivot = arr[(l + r) / 2];

  while (i <= j) {
    while (arr[i] < pivot) {
      ++i;
    }
    while (arr[j] > pivot) {
      --j;
    }
    if (i <= j) {
      std::swap(arr[i], arr[j]);
      ++i;
      --j;
    }
  }

  return {i, j};
}

void TochilinEHoarSortSimMerALL::QuickSortOMP(std::vector<int> &arr, int low, int high, int depth_limit) {
  std::vector<std::pair<int, int>> stack;
  stack.emplace_back(low, high);

  while (!stack.empty()) {
    const auto [l0, r0] = stack.back();
    stack.pop_back();

    int l = l0;
    int r = r0;

    if (l >= r) {
      continue;
    }

    const std::pair<int, int> bounds = Partition(arr, l, r);
    int i = bounds.first;
    int j = bounds.second;

    if (depth_limit > 0) {
#pragma omp task default(none) shared(arr) firstprivate(l, j, depth_limit)
      QuickSortOMP(arr, l, j, depth_limit - 1);

#pragma omp task default(none) shared(arr) firstprivate(i, r, depth_limit)
      QuickSortOMP(arr, i, r, depth_limit - 1);
    } else {
      if (l < j) {
        stack.emplace_back(l, j);
      }
      if (i < r) {
        stack.emplace_back(i, r);
      }
    }
  }

#pragma omp taskwait
}

std::vector<int> TochilinEHoarSortSimMerALL::MergeSortedVectors(const std::vector<int> &a, const std::vector<int> &b) {
  std::vector<int> result;
  result.reserve(a.size() + b.size());
  std::ranges::merge(a, b, std::back_inserter(result));
  return result;
}

bool TochilinEHoarSortSimMerALL::RunImpl() {
  auto &data = GetOutput();
  if (data.empty()) {
    return false;
  }

  int rank = 0;
  int proc_count = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &proc_count);

  const int total_size = static_cast<int>(data.size());
  const std::vector<int> counts = BuildCounts(total_size, proc_count);
  const std::vector<int> displs = BuildDisplacements(counts);

  std::vector<int> local_data(static_cast<std::size_t>(counts[static_cast<std::size_t>(rank)]), 0);
  MPI_Scatterv(data.data(), counts.data(), displs.data(), MPI_INT, local_data.data(),
               counts[static_cast<std::size_t>(rank)], MPI_INT, 0, MPI_COMM_WORLD);

  if (!local_data.empty()) {
    const int thread_count = std::max(1, ppc::util::GetNumThreads());
#pragma omp parallel default(none) shared(local_data) num_threads(thread_count) if (thread_count > 1)
    {
#pragma omp single
      QuickSortOMP(local_data, 0, static_cast<int>(local_data.size()) - 1, 3);
    }
  }

  if (rank == 0) {
    data.clear();
  } else {
    data.assign(static_cast<std::size_t>(total_size), 0);
  }

  local_data = MergeAcrossRanks(std::move(local_data), rank, proc_count);

  if (rank == 0) {
    data = std::move(local_data);
  }

  MPI_Bcast(data.data(), total_size, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

bool TochilinEHoarSortSimMerALL::PostProcessingImpl() {
  return std::ranges::is_sorted(GetOutput());
}

}  // namespace tochilin_e_hoar_sort_sim_mer
