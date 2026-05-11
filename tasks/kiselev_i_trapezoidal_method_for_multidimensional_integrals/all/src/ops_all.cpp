#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/all/include/ops_all.hpp"

#include <mpi.h>

#include <cmath>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/common/include/common.hpp"

namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals {

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

double KiselevITestTaskALL::FunctionTypeChoose(int type_x, double x, double y) {
  switch (type_x) {
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

double KiselevITestTaskALL::ComputeIntegral(const std::vector<int> &steps) {
  const auto &in = GetInput();

  const int nx = steps[0];
  const int ny = steps[1];

  const double lx = in.left_bounds[0];
  const double ly = in.left_bounds[1];

  const double rx = in.right_bounds[0];
  const double ry = in.right_bounds[1];

  const double hx = (rx - lx) / static_cast<double>(nx);

  const double hy = (ry - ly) / static_cast<double>(ny);

  int mpi_initialized = 0;

  MPI_Initialized(&mpi_initialized);

  int rank = 0;
  int size = 1;

  if (mpi_initialized) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }

  const int total_rows = nx + 1;

  const int rows_per_proc = total_rows / size;

  const int remainder = total_rows % size;

  const int begin = rank * rows_per_proc + std::min(rank, remainder);

  const int local_rows = rows_per_proc + (rank < remainder ? 1 : 0);

  const int end = begin + local_rows;

  double local_sum = 0.0;

  switch (in.type_function) {
    case 0: {
#pragma omp parallel for collapse(2) reduction(+ : local_sum) schedule(static)
      for (int i = begin; i < end; i++) {
        for (int j = 0; j <= ny; j++) {
          const double x = lx + (i * hx);
          const double y = ly + (j * hy);

          const double wx = (i == 0 || i == nx) ? 0.5 : 1.0;

          const double wy = (j == 0 || j == ny) ? 0.5 : 1.0;

          local_sum += wx * wy * ((x * x) + (y * y));
        }
      }

      break;
    }

    case 1: {
#pragma omp parallel for collapse(2) reduction(+ : local_sum) schedule(static)
      for (int i = begin; i < end; i++) {
        for (int j = 0; j <= ny; j++) {
          const double x = lx + (i * hx);
          const double y = ly + (j * hy);

          const double wx = (i == 0 || i == nx) ? 0.5 : 1.0;

          const double wy = (j == 0 || j == ny) ? 0.5 : 1.0;

          local_sum += wx * wy * (std::sin(x) * std::cos(y));
        }
      }

      break;
    }

    case 2: {
#pragma omp parallel for collapse(2) reduction(+ : local_sum) schedule(static)
      for (int i = begin; i < end; i++) {
        for (int j = 0; j <= ny; j++) {
          const double x = lx + (i * hx);
          const double y = ly + (j * hy);

          const double wx = (i == 0 || i == nx) ? 0.5 : 1.0;

          const double wy = (j == 0 || j == ny) ? 0.5 : 1.0;

          local_sum += wx * wy * (std::sin(x) + std::cos(y));
        }
      }

      break;
    }

    case 3: {
#pragma omp parallel for collapse(2) reduction(+ : local_sum) schedule(static)
      for (int i = begin; i < end; i++) {
        for (int j = 0; j <= ny; j++) {
          const double x = lx + (i * hx);
          const double y = ly + (j * hy);

          const double wx = (i == 0 || i == nx) ? 0.5 : 1.0;

          const double wy = (j == 0 || j == ny) ? 0.5 : 1.0;

          local_sum += wx * wy * std::exp(x + y);
        }
      }

      break;
    }

    default: {
#pragma omp parallel for collapse(2) reduction(+ : local_sum) schedule(static)
      for (int i = begin; i < end; i++) {
        for (int j = 0; j <= ny; j++) {
          const double x = lx + (i * hx);
          const double y = ly + (j * hy);

          const double wx = (i == 0 || i == nx) ? 0.5 : 1.0;

          const double wy = (j == 0 || j == ny) ? 0.5 : 1.0;

          local_sum += wx * wy * (x + y);
        }
      }

      break;
    }
  }

  double global_sum = local_sum;

  if (mpi_initialized && size > 1) {
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }

  return global_sum * hx * hy;
}

bool KiselevITestTaskALL::RunImpl() {
  std::vector<int> steps = GetInput().step_n_size;

  const auto &in = GetInput();

  if (in.left_bounds.size() != 2 || in.right_bounds.size() != 2 || in.step_n_size.size() != 2) {
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

  int iter = 0;

  const int max_iter = 1;

  while (iter < max_iter) {
    for (auto &s : steps) {
      s *= 2;
    }

    current = ComputeIntegral(steps);

    if (std::abs(current - prev) < epsilon) {
      break;
    }

    prev = current;

    iter++;
  }

  GetOutput() = current;

  return true;
}

bool KiselevITestTaskALL::PostProcessingImpl() {
  return true;
}

}  // namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals
