#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/common/include/common.hpp"

namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals {

namespace {

double EvaluateFunction(int type_function, double x, double y) {
  switch (type_function) {
    case 0:
      return (x * x) + (y * y);

    case 1:
      return std::sin(x) * std::cos(y);

    case 2:
      return std::sin(x) + std::cos(y);

    case 3:
      return std::exp(x + y);

    default:
      return x + y;
  }
}

double ComputeLocalIntegral(const InType &in, const std::vector<int> &steps, int begin, int end) {
  const int nx = steps[0];
  const int ny = steps[1];

  const double hx = (in.right_bounds[0] - in.left_bounds[0]) / static_cast<double>(nx);
  const double hy = (in.right_bounds[1] - in.left_bounds[1]) / static_cast<double>(ny);

  double local_sum = 0.0;

#pragma omp parallel for reduction(+ : local_sum) schedule(static) default(none) \
    shared(in, steps, begin, end, hx, hy, nx, ny)
  for (int i = begin; i < end; ++i) {
    for (int j = 0; j <= ny; ++j) {
      const double x = in.left_bounds[0] + (static_cast<double>(i) * hx);
      const double y = in.left_bounds[1] + (static_cast<double>(j) * hy);
      const double wx = ((i == 0) || (i == nx)) ? 0.5 : 1.0;
      const double wy = ((j == 0) || (j == ny)) ? 0.5 : 1.0;

      local_sum += wx * wy * EvaluateFunction(in.type_function, x, y);
    }
  }

  return local_sum * hx * hy;
}

}  // namespace

KiselevITestTaskALL::KiselevITestTaskALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());

  GetInput() = in;

  GetOutput() = 0.0;
}

bool KiselevITestTaskALL::ValidationImpl() {
  return true;
}

bool KiselevITestTaskALL::PreProcessingImpl() {
  GetOutput() = 0.0;

  return true;
}

double KiselevITestTaskALL::ComputeIntegral(const std::vector<int> &steps) {
  const auto &in = GetInput();

  int mpi_initialized = 0;

  MPI_Initialized(&mpi_initialized);

  int rank = 0;
  int size = 1;

  if (mpi_initialized != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }

  const int total_rows = steps[0] + 1;
  const int rows_per_proc = total_rows / size;
  const int remainder = total_rows % size;
  const int begin = (rank * rows_per_proc) + std::min(rank, remainder);
  const int local_rows = rows_per_proc + ((rank < remainder) ? 1 : 0);
  const int end = begin + local_rows;

  double local_result = ComputeLocalIntegral(in, steps, begin, end);
  double global_result = local_result;

  if ((mpi_initialized != 0) && (size > 1)) {
    MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }

  return global_result;
}

bool KiselevITestTaskALL::RunImpl() {
  std::vector<int> steps = GetInput().step_n_size;

  const auto &in = GetInput();

  if ((in.left_bounds.size() != 2) || (in.right_bounds.size() != 2) || (in.step_n_size.size() != 2)) {
    GetOutput() = 0.0;

    return true;
  }

  const double epsilon = in.epsilon;

  if (epsilon <= 0.0) {
    GetOutput() = ComputeIntegral(steps);

    return true;
  }

  double prev = ComputeIntegral(steps);

  double current = prev;

  constexpr int kMaxIter = 1;

  for (int iter = 0; iter < kMaxIter; ++iter) {
    for (auto &s : steps) {
      s *= 2;
    }

    current = ComputeIntegral(steps);

    if (std::abs(current - prev) < epsilon) {
      break;
    }

    prev = current;
  }

  GetOutput() = current;

  return true;
}

bool KiselevITestTaskALL::PostProcessingImpl() {
  return true;
}

}  // namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals
