#include "trofimov_n_hoar_sort_batcher/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <utility>
#include <vector>

#include "oneapi/tbb/parallel_sort.h"
#include "trofimov_n_hoar_sort_batcher/common/include/common.hpp"
#include "util/include/util.hpp"

namespace trofimov_n_hoar_sort_batcher {

namespace {

constexpr int kQuickSortDepthLimit = 4;
constexpr int kSequentialThreshold = 1024;

struct MpiDistribution {
  int rank = 0;
  int size = 1;
  int n = 0;
  std::vector<int> counts;
  std::vector<int> displs;
};

int HoarePartition(std::vector<int> &arr, int left, int right) {
  const int pivot = arr[left + ((right - left) / 2)];
  int i = left - 1;
  int j = right + 1;

  while (true) {
    while (arr[++i] < pivot) {
    }
    while (arr[--j] > pivot) {
    }
    if (i >= j) {
      return j;
    }
    std::swap(arr[i], arr[j]);
  }
}

void TbbSortRange(std::vector<int> &arr, int left, int right) {
  if (left >= right) {
    return;
  }
  oneapi::tbb::parallel_sort(arr.begin() + left, arr.begin() + right + 1);
}

void QuickSortOmpTask(std::vector<int> &arr, int left, int right, int depth_limit) {
  if (left >= right) {
    return;
  }

  if ((right - left) < kSequentialThreshold || depth_limit <= 0) {
    TbbSortRange(arr, left, right);
    return;
  }

  const int split = HoarePartition(arr, left, right);

#pragma omp task default(none) shared(arr) \
    firstprivate(left, split, depth_limit) if ((split - left) > kSequentialThreshold)
  QuickSortOmpTask(arr, left, split, depth_limit - 1);

#pragma omp task default(none) shared(arr) \
    firstprivate(split, right, depth_limit) if ((right - (split + 1)) > kSequentialThreshold)
  QuickSortOmpTask(arr, split + 1, right, depth_limit - 1);

#pragma omp taskwait
}

MpiDistribution BuildDistribution(int n) {
  MpiDistribution dist;
  dist.n = n;
  MPI_Comm_rank(MPI_COMM_WORLD, &dist.rank);
  MPI_Comm_size(MPI_COMM_WORLD, &dist.size);

  dist.counts.assign(dist.size, 0);
  dist.displs.assign(dist.size, 0);

  const int base = n / dist.size;
  const int rem = n % dist.size;

  for (int i = 0; i < dist.size; i++) {
    dist.counts[i] = base + (i < rem ? 1 : 0);
    if (i > 0) {
      dist.displs[i] = dist.displs[i - 1] + dist.counts[i - 1];
    }
  }

  return dist;
}

std::vector<int> ScatterData(const std::vector<int> &data, const MpiDistribution &dist) {
  std::vector<int> local(dist.counts[dist.rank]);
  MPI_Scatterv(data.data(), dist.counts.data(), dist.displs.data(), MPI_INT, local.data(), dist.counts[dist.rank],
               MPI_INT, 0, MPI_COMM_WORLD);
  return local;
}

void SortLocalChunkHybrid(std::vector<int> &local) {
  if (local.empty()) {
    return;
  }

#pragma omp parallel num_threads(ppc::util::GetNumThreads()) default(none) shared(local)
  {
#pragma omp single nowait
    {
      QuickSortOmpTask(local, 0, static_cast<int>(local.size()) - 1, kQuickSortDepthLimit);
    }
  }
}

std::vector<int> GatherChunks(const std::vector<int> &local, const MpiDistribution &dist) {
  std::vector<int> gathered;
  if (dist.rank == 0) {
    gathered.resize(dist.n);
  }

  MPI_Gatherv(local.data(), dist.counts[dist.rank], MPI_INT, gathered.data(), dist.counts.data(), dist.displs.data(),
              MPI_INT, 0, MPI_COMM_WORLD);

  return gathered;
}

void MergeSortedChunksOnRoot(std::vector<int> &gathered, const MpiDistribution &dist) {
  if (dist.rank != 0 || gathered.empty()) {
    return;
  }

  int merged_size = dist.counts[0];
  for (int proc = 1; proc < dist.size; proc++) {
    const int right_end = merged_size + dist.counts[proc];
    std::inplace_merge(gathered.begin(), gathered.begin() + merged_size, gathered.begin() + right_end);
    merged_size = right_end;
  }
}

void BroadcastResult(std::vector<int> &data, const MpiDistribution &dist) {
  MPI_Bcast(data.data(), dist.n, MPI_INT, 0, MPI_COMM_WORLD);
}

}  // namespace

TrofimovNHoarSortBatcherALL::TrofimovNHoarSortBatcherALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TrofimovNHoarSortBatcherALL::ValidationImpl() {
  return true;
}

bool TrofimovNHoarSortBatcherALL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool TrofimovNHoarSortBatcherALL::RunImpl() {
  auto &data = GetOutput();
  if (data.size() <= 1) {
    return true;
  }

  const auto dist = BuildDistribution(static_cast<int>(data.size()));
  auto local = ScatterData(data, dist);
  SortLocalChunkHybrid(local);

  auto gathered = GatherChunks(local, dist);
  MergeSortedChunksOnRoot(gathered, dist);

  if (dist.rank == 0) {
    data = std::move(gathered);
  }

  BroadcastResult(data, dist);
  return true;
}

bool TrofimovNHoarSortBatcherALL::PostProcessingImpl() {
  return true;
}

}  // namespace trofimov_n_hoar_sort_batcher
