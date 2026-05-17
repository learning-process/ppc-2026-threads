#include "tsibareva_e_integral_calculate_trapezoid_method/all/include/ops_all.hpp"

#include <mpi.h>

#include <cmath>
#include <functional>
#include <vector>

#include "tsibareva_e_integral_calculate_trapezoid_method/common/include/common.hpp"

namespace tsibareva_e_integral_calculate_trapezoid_method {

TsibarevaEIntegralCalculateTrapezoidMethodALL::TsibarevaEIntegralCalculateTrapezoidMethodALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodALL::ValidationImpl() {
  return true;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodALL::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

double TsibarevaEIntegralCalculateTrapezoidMethodALL::ComputePartialSum(
    int begin, int finish, const std::vector<double> &lo, const std::vector<double> &h, const std::vector<int> &sizes,
    const std::vector<int> &steps, int dim, const std::function<double(const std::vector<double> &)> &f) {
  double partial = 0.0;
#pragma omp parallel default(none) shared(begin, finish, dim, h, sizes, lo, steps, f) reduction(+ : partial)
  {
    std::vector<double> point(dim);
#pragma omp for
    for (int node = begin; node < finish; ++node) {
      int remainder_idx = node;
      double node_weight = 1.0;
      for (int i = dim - 1; i >= 0; --i) {
        int idx = remainder_idx % sizes[i];
        remainder_idx /= sizes[i];
        if (idx == 0 || idx == steps[i]) {
          node_weight *= 0.5;
        }
        point[i] = lo[i] + (idx * h[i]);
      }
      partial += node_weight * f(point);
    }
  }
  return partial;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodALL::RunImpl() {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const auto &lo = GetInput().lo;
  const auto &hi = GetInput().hi;
  const auto &steps = GetInput().steps;
  const auto &f = GetInput().f;
  int dim = GetInput().dim;

  MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
  std::vector<double> lo_vec(dim);
  std::vector<double> hi_vec(dim);
  std::vector<int> steps_vec(dim);
  if (rank == 0) {
    lo_vec = lo;
    hi_vec = hi;
    steps_vec = steps;
  }
  MPI_Bcast(lo_vec.data(), dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(hi_vec.data(), dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(steps_vec.data(), dim, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<double> h(dim);
  std::vector<int> sizes(dim);
  int total_nodes = 1;
  for (int i = 0; i < dim; ++i) {
    h[i] = (hi_vec[i] - lo_vec[i]) / static_cast<double>(steps_vec[i]);
    sizes[i] = steps_vec[i] + 1;
    total_nodes *= sizes[i];
  }

  std::vector<int> all_starts(size);
  std::vector<int> all_ends(size);
  if (rank == 0) {
    int nodes_per_proc = total_nodes / size;
    int remainder = total_nodes % size;
    int start = 0;
    for (int i = 0; i < size; ++i) {
      all_starts[i] = start;
      int end = start + nodes_per_proc + (i < remainder ? 1 : 0);
      all_ends[i] = end;
      start = end;
    }
  }

  int my_start = 0;
  int my_end = 0;
  MPI_Scatter(all_starts.data(), 1, MPI_INT, &my_start, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(all_ends.data(), 1, MPI_INT, &my_end, 1, MPI_INT, 0, MPI_COMM_WORLD);

  double local_sum = ComputePartialSum(my_start, my_end, lo_vec, h, sizes, steps_vec, dim, f);

  double global_sum = 0.0;
  MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    double res_h = 1.0;
    for (int i = 0; i < dim; ++i) {
      res_h *= h[i];
    }
    GetOutput() = global_sum * res_h;
  }

  return true;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodALL::PostProcessingImpl() {
  MPI_Bcast(&GetOutput(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

}  // namespace tsibareva_e_integral_calculate_trapezoid_method
