#include "../include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <tuple>
#include <utility>
#include <vector>

#include "../../common/include/common.hpp"
#include "util/include/util.hpp"

namespace kutergin_a_multidim_trapezoid {

KuterginAMultidimTrapezoidALL::KuterginAMultidimTrapezoidALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool KuterginAMultidimTrapezoidALL::ValidationImpl() {
  const auto &input = GetInput();
  const auto &limits = std::get<1>(input);
  const int n_steps = std::get<2>(input);

  if (limits.empty()) {
    return false;
  }
  return n_steps > 0;
}

bool KuterginAMultidimTrapezoidALL::PreProcessingImpl() {
  local_input_ = GetInput();
  res_ = 0.0;
  return true;
}

bool KuterginAMultidimTrapezoidALL::RunImpl() {
  int rank = 0;
  int size = 1;
  const bool is_mpi = ppc::util::IsUnderMpirun();

  if (is_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }

  size_t dims = (rank == 0) ? std::get<1>(local_input_).size() : 0;

  DistributeData(rank, dims);

  const auto &func = std::get<0>(local_input_);
  const auto &limits = std::get<1>(local_input_);
  const int n_steps = std::get<2>(local_input_);

  size_t total_nodes = 1;
  std::vector<double> h(dims);
  double cell_volume = 1.0;

  for (size_t i = 0; i < dims; ++i) {
    total_nodes *= (static_cast<size_t>(n_steps) + 1);
    h[i] = (limits[i].second - limits[i].first) / n_steps;
    cell_volume *= h[i];
  }

  const size_t proc_chunk = total_nodes / size;
  const size_t proc_remainder = total_nodes % size;

  const size_t my_proc_count = proc_chunk + (std::cmp_less(rank, proc_remainder) ? 1 : 0);
  const size_t my_proc_start = (rank * proc_chunk) + std::min(static_cast<size_t>(rank), proc_remainder);

  double local_sum = 0.0;

  if (my_proc_count > 0) {
    int num_threads = ppc::util::GetNumThreads();
    omp_set_num_threads(num_threads);

#pragma omp parallel default(none) shared(h, dims, my_proc_start, my_proc_count, func, n_steps, limits) \
    reduction(+ : local_sum)
    {
      int tid = omp_get_thread_num();
      int t_count = omp_get_num_threads();

      size_t thread_chunk = my_proc_count / t_count;
      size_t thread_remainder = my_proc_count % t_count;

      size_t my_thread_count = thread_chunk + (std::cmp_less(tid, thread_remainder) ? 1 : 0);
      size_t my_thread_start =
          my_proc_start + (tid * thread_chunk) + std::min(static_cast<size_t>(tid), thread_remainder);

      local_sum += CalculateChunkSum(my_thread_start, my_thread_start + my_thread_count, h, limits, n_steps, func);
    }
  }

  if (is_mpi) {
    double global_sum = 0.0;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
      res_ = global_sum * cell_volume;
    }
  } else {
    res_ = local_sum * cell_volume;
  }

  return true;
}

bool KuterginAMultidimTrapezoidALL::PostProcessingImpl() {
  if (ppc::util::IsUnderMpirun()) {
    MPI_Bcast(&res_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  GetOutput() = res_;
  return true;
}

void KuterginAMultidimTrapezoidALL::DistributeData(int rank, size_t &dims) {
  if (!ppc::util::IsUnderMpirun()) {
    return;
  }

  int dims_io = static_cast<int>(dims);
  MPI_Bcast(&dims_io, 1, MPI_INT, 0, MPI_COMM_WORLD);
  dims = static_cast<size_t>(dims_io);

  auto &limits = std::get<1>(local_input_);
  auto &n_steps = std::get<2>(local_input_);

  if (rank != 0) {
    limits.resize(dims);
  }

  for (size_t i = 0; i < dims; ++i) {
    MPI_Bcast(&limits[i].first, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&limits[i].second, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  MPI_Bcast(&n_steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

double KuterginAMultidimTrapezoidALL::CalculateChunkSum(
    size_t start_idx, size_t end_idx, const std::vector<double> &h,
    const std::vector<std::pair<double, double>> &limits, int n_steps,
    const std::function<double(const std::vector<double> &)> &func) {
  if (start_idx >= end_idx) {
    return 0.0;
  }

  const size_t dims = limits.size();
  std::vector<double> coords(dims);
  double chunk_sum = 0.0;

  for (size_t i = start_idx; i < end_idx; ++i) {
    double weight = 1.0;
    size_t temp_idx = i;

    for (size_t dim = 0; dim < dims; ++dim) {
      const int nodes_in_dim = n_steps + 1;
      const int coord_idx = static_cast<int>(temp_idx % nodes_in_dim);
      temp_idx /= nodes_in_dim;

      coords[dim] = limits[dim].first + (static_cast<double>(coord_idx) * h[dim]);

      if (coord_idx == 0 || coord_idx == n_steps) {
        weight *= 0.5;
      }
    }
    chunk_sum += weight * func(coords);
  }

  return chunk_sum;
}

}  // namespace kutergin_a_multidim_trapezoid
