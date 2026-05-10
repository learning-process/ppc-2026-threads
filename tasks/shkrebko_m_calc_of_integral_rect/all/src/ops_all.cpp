#include "shkrebko_m_calc_of_integral_rect/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
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
  for (size_t i = 0; i < input.n_steps.size(); ++i) {
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

  size_t dims = local_input_.limits.size();
  int dims_int = static_cast<int>(dims);
  MPI_Bcast(&dims_int, 1, MPI_INT, 0, MPI_COMM_WORLD);
  dims = static_cast<size_t>(dims_int);

  if (rank != 0) {
    local_input_.limits.resize(dims);
    local_input_.n_steps.resize(dims);
  }

  for (size_t i = 1; i < dims; ++i) {
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

  int base = global_steps / size;
  int remainder = global_steps % size;
  local_steps = base + (rank < remainder ? 1 : 0);
  local_offset = rank * base + std::min(rank, remainder);

  double step = (global_right - global_left) / global_steps;
  local_left = global_left + local_offset * step;
  local_right = local_left + local_steps * step;
}

double ShkrebkoMCalcOfIntegralRectALL::ComputeSliceSum(double left0, double right0, int steps0,
                                                       const std::vector<double> &h_other,
                                                       const std::vector<std::pair<double, double>> &limits_other,
                                                       const std::vector<int> &n_steps_other) const {
  const size_t other_dims = limits_other.size();
  // Если других измерений нет, то просто интегрируем по x0
  if (other_dims == 0) {
    double h0 = (right0 - left0) / steps0;
    double sum = 0.0;
    for (int i = 0; i < steps0; ++i) {
      double x0 = left0 + (i + 0.5) * h0;
      sum += local_input_.func({x0});
    }
    return sum * h0;
  }

  double h0 = (right0 - left0) / steps0;
  double total = 0.0;

#pragma omp parallel for reduction(+ : total) schedule(static)
  for (int i = 0; i < steps0; ++i) {
    double x0 = left0 + (i + 0.5) * h0;

    std::vector<int> indices(other_dims, 0);
    double local_sum = 0.0;
    size_t total_other_points = 1;
    for (size_t d = 0; d < other_dims; ++d) {
      total_other_points *= n_steps_other[d];
    }

    for (size_t idx = 0; idx < total_other_points; ++idx) {
      std::vector<double> point(other_dims + 1);
      point[0] = x0;
      size_t tmp = idx;
      for (int d = static_cast<int>(other_dims) - 1; d >= 0; --d) {
        int coord = static_cast<int>(tmp % n_steps_other[d]);
        tmp /= n_steps_other[d];
        point[d + 1] = limits_other[d].first + (coord + 0.5) * h_other[d];
      }
      local_sum += local_input_.func(point);
    }
    total += local_sum;
  }

  double vol_other = 1.0;
  for (double h : h_other) {
    vol_other *= h;
  }
  return total * h0 * vol_other;
}

bool ShkrebkoMCalcOfIntegralRectALL::RunImpl() {
  int rank = 0, size = 1;
  bool is_mpi = ppc::util::IsUnderMpirun();
  if (is_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }

  BroadcastCommonData(rank);

  const size_t dims = local_input_.limits.size();
  if (dims == 0) {
    return false;
  }
  if (!local_input_.func) {
    return false;
  }

  double my_left, my_right;
  int my_steps, my_offset;
  AssignMpiSlice(rank, size, my_left, my_right, my_steps, my_offset);

  std::vector<double> h_other;
  std::vector<std::pair<double, double>> limits_other;
  std::vector<int> n_steps_other;
  for (size_t i = 1; i < dims; ++i) {
    limits_other.push_back(local_input_.limits[i]);
    n_steps_other.push_back(local_input_.n_steps[i]);
    double h = (local_input_.limits[i].second - local_input_.limits[i].first) / local_input_.n_steps[i];
    h_other.push_back(h);
  }

  double local_slice_integral = ComputeSliceSum(my_left, my_right, my_steps, h_other, limits_other, n_steps_other);

  double global_integral = 0.0;
  if (is_mpi) {
    MPI_Reduce(&local_slice_integral, &global_integral, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
      res_ = global_integral;
    }
  } else {
    res_ = local_slice_integral;
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
