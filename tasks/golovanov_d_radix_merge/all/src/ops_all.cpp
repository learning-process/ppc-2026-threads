#include "golovanov_d_radix_merge/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <utility>
#include <vector>

#include "golovanov_d_radix_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace {

inline std::uint64_t DoubleToSortableUint64(double val) {
  std::uint64_t bits = 0;
  std::memcpy(&bits, &val, sizeof(double));
  const std::uint64_t sign_mask = std::uint64_t{1} << 63;
  return ((bits & sign_mask) != 0U) ? ~bits : bits ^ sign_mask;
}

void RadixSortStable(std::vector<double> &data, int /*num_threads*/) {
  if (data.size() < 2) {
    return;
  }
  std::ranges::stable_sort(data,
                           [](double a, double b) { return DoubleToSortableUint64(a) < DoubleToSortableUint64(b); });
}

void MergeSortedVectors(std::vector<double> &left, const std::vector<double> &right) {
  if (right.empty()) {
    return;
  }
  if (left.empty()) {
    left = right;
    return;
  }
  std::vector<double> merged(left.size() + right.size());
  std::ranges::merge(left, right, merged.begin(),
                     [](double a, double b) { return DoubleToSortableUint64(a) < DoubleToSortableUint64(b); });
  left.swap(merged);
}

struct DistributionParams {
  std::vector<int> send_counts;
  std::vector<int> displs;
  int global_size_int{};
};

DistributionParams ComputeDistribution(int world_size, uint64_t global_size) {
  DistributionParams params;
  params.global_size_int = static_cast<int>(global_size);
  params.send_counts.resize(world_size, 0);
  params.displs.resize(world_size, 0);

  const int base_count = params.global_size_int / world_size;
  const int remainder = params.global_size_int % world_size;
  int offset = 0;

  for (int proc = 0; proc < world_size; ++proc) {
    params.send_counts[proc] = base_count + (proc < remainder ? 1 : 0);
    params.displs[proc] = offset;
    offset += params.send_counts[proc];
  }

  return params;
}

bool PerformMergeStep(std::vector<double> &local_data, int rank, int step, int world_size) {
  const int group_size = step << 1;

  if ((rank % group_size) == 0) {
    const int src = rank + step;
    if (src < world_size) {
      int recv_count = 0;
      MPI_Recv(&recv_count, 1, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      if (recv_count > 0) {
        std::vector<double> received_data(recv_count);
        MPI_Recv(received_data.data(), recv_count, MPI_DOUBLE, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MergeSortedVectors(local_data, received_data);
      }
    }
    return true;
  }

  const int dst = rank - step;
  const int send_count = static_cast<int>(local_data.size());
  MPI_Send(&send_count, 1, MPI_INT, dst, 0, MPI_COMM_WORLD);

  if (send_count > 0) {
    MPI_Send(local_data.data(), send_count, MPI_DOUBLE, dst, 1, MPI_COMM_WORLD);
  }
  local_data.clear();
  return false;
}

void ParallelMergeSort(std::vector<double> &local_data, int rank, int world_size, int num_threads) {
  RadixSortStable(local_data, num_threads);

  for (int step = 1; step < world_size; step <<= 1) {
    if (!PerformMergeStep(local_data, rank, step, world_size)) {
      break;
    }
  }
}

void GatherResults(std::vector<double> &local_data, int rank, int global_size, std::vector<double> &result) {
  result.resize(global_size);
  if (rank == 0) {
    result.swap(local_data);
  }
  MPI_Bcast(result.data(), global_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

}  // namespace

namespace golovanov_d_radix_merge {

GolovanovDRadixMergeALL::GolovanovDRadixMergeALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool GolovanovDRadixMergeALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int is_valid = 1;
  if (rank == 0) {
    is_valid = !GetInput().empty() ? 1 : 0;
  }
  MPI_Bcast(&is_valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return is_valid != 0;
}

bool GolovanovDRadixMergeALL::PreProcessingImpl() {
  return true;
}

bool GolovanovDRadixMergeALL::RunImpl() {
  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  const int num_threads = ppc::util::GetNumThreads();

  uint64_t global_size = 0;
  if (rank == 0) {
    global_size = GetInput().size();
  }
  MPI_Bcast(&global_size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  DistributionParams dist_params = ComputeDistribution(world_size, global_size);

  std::vector<double> local_data(dist_params.send_counts[rank]);
  double *send_buffer = (rank == 0) ? GetInput().data() : nullptr;

  MPI_Scatterv(send_buffer, dist_params.send_counts.data(), dist_params.displs.data(), MPI_DOUBLE,
               local_data.empty() ? nullptr : local_data.data(), dist_params.send_counts[rank], MPI_DOUBLE, 0,
               MPI_COMM_WORLD);

  ParallelMergeSort(local_data, rank, world_size, num_threads);

  std::vector<double> result;
  GatherResults(local_data, rank, dist_params.global_size_int, result);

  GetOutput() = std::move(result);
  return true;
}

bool GolovanovDRadixMergeALL::PostProcessingImpl() {
  return true;
}

}  // namespace golovanov_d_radix_merge
