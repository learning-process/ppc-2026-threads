#include "shkrebko_m_calc_of_integral_rect/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "shkrebko_m_calc_of_integral_rect/common/include/common.hpp"
#include "util/include/util.hpp"

namespace shkrebko_m_calc_of_integral_rect {

ShkrebkoMCalcOfIntegralRectALL::ShkrebkoMCalcOfIntegralRectALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ShkrebkoMCalcOfIntegralRectALL::ValidationImpl() {
  const auto &input = GetInput();
  if (!input.func || input.limits.empty() || input.limits.size() != input.n_steps.size()) {
    return false;
  }
  for (std::size_t i = 0; i < input.n_steps.size(); ++i) {
    if (input.n_steps[i] <= 0 || input.limits[i].first >= input.limits[i].second) {
      return false;
    }
  }
  return true;
}

bool ShkrebkoMCalcOfIntegralRectALL::PreProcessingImpl() {
  local_input_ = GetInput();
  res_ = 0.0;
  return true;
}

void ShkrebkoMCalcOfIntegralRectALL::BroadcastCommonData(int rank) {
  if (!ppc::util::IsUnderMpirun()) {
    return;
  }

  std::size_t dims = local_input_.limits.size();
  int dims_int = static_cast<int>(dims);
  MPI_Bcast(&dims_int, 1, MPI_INT, 0, MPI_COMM_WORLD);
  dims = static_cast<std::size_t>(dims_int);

  if (rank != 0) {
    local_input_.limits.resize(dims);
    local_input_.n_steps.resize(dims);
  }

  for (std::size_t i = 1; i < dims; ++i) {
    MPI_Bcast(&local_input_.limits[i].first, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&local_input_.limits[i].second, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&local_input_.n_steps[i], 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
}

void ShkrebkoMCalcOfIntegralRectALL::AssignMpiSlice(int rank, int size, double &local_left, double &local_right,
                                                    int &local_steps, int &local_offset) {
  const double global_left = local_input_.limits[0].first;
  const double global_right = local_input_.limits[0].second;
  const int global_steps = local_input_.n_steps[0];

  const int base = global_steps / size;
  const int remainder = global_steps % size;
  local_steps = base + (rank < remainder ? 1 : 0);
  local_offset = (rank * base) + std::min(rank, remainder);

  const double step = (global_right - global_left) / static_cast<double>(global_steps);
  local_left = global_left + (static_cast<double>(local_offset) * step);
  local_right = local_left + (static_cast<double>(local_steps) * step);
}

double ShkrebkoMCalcOfIntegralRectALL::ComputeSliceSum(double left0, double right0, int steps0,
                                                       const std::vector<double> &h_other,
                                                       const std::vector<std::pair<double, double>> &limits_other,
                                                       const std::vector<int> &n_steps_other) const {
  const std::size_t other_dims = limits_other.size();

  if (other_dims == 0) {
    const double h0 = (right0 - left0) / static_cast<double>(steps0);
    double sum = 0.0;
    for (int i = 0; i < steps0; ++i) {
      const double x0 = left0 + ((static_cast<double>(i) + 0.5) * h0);
      sum += local_input_.func({x0});
    }
    return sum * h0;
  }

  const double h0 = (right0 - left0) / static_cast<double>(steps0);
  double total = 0.0;

#pragma omp parallel for default(none) shared(left0, h0, steps0, h_other, limits_other, n_steps_other, other_dims) \
    reduction(+ : total) schedule(static)
  for (int i = 0; i < steps0; ++i) {
    const double x0 = left0 + ((static_cast<double>(i) + 0.5) * h0);

    std::vector<int> indices(other_dims, 0);
    double local_sum = 0.0;

    std::size_t total_other_points = 1;
    for (std::size_t dim = 0; dim < other_dims; ++dim) {
      total_other_points *= static_cast<std::size_t>(n_steps_other[dim]);
    }

    for (std::size_t idx = 0; idx < total_other_points; ++idx) {
      std::vector<double> point(other_dims + 1);
      point[0] = x0;
      std::size_t tmp = idx;
      for (int dim = static_cast<int>(other_dims) - 1; dim >= 0; --dim) {
        const int coord = static_cast<int>(tmp % static_cast<std::size_t>(n_steps_other[dim]));
        tmp /= static_cast<std::size_t>(n_steps_other[dim]);
        point[dim + 1] = limits_other[dim].first + ((static_cast<double>(coord) + 0.5) * h_other[dim]);
      }
      local_sum += local_input_.func(point);
    }
    total += local_sum;
  }

  double volume_other = 1.0;
  for (double h : h_other) {
    volume_other *= h;
  }
  return total * h0 * volume_other;
}

bool ShkrebkoMCalcOfIntegralRectALL::RunImpl() {
  int rank = 0;
  int size = 1;
  const bool is_mpi = ppc::util::IsUnderMpirun();
  if (is_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }

  BroadcastCommonData(rank);

  const std::size_t dims = local_input_.limits.size();
  if (dims == 0) {
    return false;
  }
  if (!local_input_.func) {
    return false;
  }

  double my_left = 0.0;
  double my_right = 0.0;
  int my_steps = 0;
  int my_offset = 0;
  AssignMpiSlice(rank, size, my_left, my_right, my_steps, my_offset);

  std::vector<double> h_other;
  std::vector<std::pair<double, double>> limits_other;
  std::vector<int> n_steps_other;
  for (std::size_t i = 1; i < dims; ++i) {
    limits_other.push_back(local_input_.limits[i]);
    n_steps_other.push_back(local_input_.n_steps[i]);
    const double h =
        (local_input_.limits[i].second - local_input_.limits[i].first) / static_cast<double>(local_input_.n_steps[i]);
    h_other.push_back(h);
  }

  const double local_slice = ComputeSliceSum(my_left, my_right, my_steps, h_other, limits_other, n_steps_other);

  double global_integral = 0.0;
  if (is_mpi) {
    MPI_Reduce(&local_slice, &global_integral, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
      res_ = global_integral;
    }
  } else {
    res_ = local_slice;
  }

  return true;
}

bool ShkrebkoMCalcOfIntegralRectALL::PostProcessingImpl() {
  int rank = 0;
  if (ppc::util::IsUnderMpirun()) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }
  if (rank == 0) {
    GetOutput() = res_;
  }
  return true;
}

}  // namespace shkrebko_m_calc_of_integral_rect
