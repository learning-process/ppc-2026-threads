#include "dergynov_s_integrals_multistep_rectangle/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>

#include "dergynov_s_integrals_multistep_rectangle/common/include/common.hpp"
#include "util/include/util.hpp"

namespace dergynov_s_integrals_multistep_rectangle {

DergynovSIntegralsMultistepRectangleALL::DergynovSIntegralsMultistepRectangleALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool DergynovSIntegralsMultistepRectangleALL::ValidationImpl() {
  const auto &[func, borders, n] = GetInput();

  if (borders.empty()) {
    return false;
  }
  for (const auto &[left, right] : borders) {
    if (!std::isfinite(left) || !std::isfinite(right)) {
      return false;
    }
    if (left >= right) {
      return false;
    }
  }

  return func && (n > 0);
}

bool DergynovSIntegralsMultistepRectangleALL::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

namespace {

void FillPoint(size_t linear_idx, int n, int dim, const std::vector<std::pair<double, double>> &borders,
               const std::vector<double> &h, std::vector<double> &point) {
  size_t tmp = linear_idx;
  for (int axis = dim - 1; axis >= 0; --axis) {
    int idx_val = static_cast<int>(tmp % static_cast<size_t>(n));
    tmp /= static_cast<size_t>(n);
    point[axis] = borders[axis].first + ((static_cast<double>(idx_val) + 0.5) * h[axis]);
  }
}

}  // namespace

bool DergynovSIntegralsMultistepRectangleALL::RunImpl() {
  const auto &input = GetInput();
  const auto &func = std::get<0>(input);
  const auto &borders = std::get<1>(input);
  const int n = std::get<2>(input);
  const int dim = static_cast<int>(borders.size());

  std::vector<double> h(dim);
  double cell_volume = 1.0;
  for (int i = 0; i < dim; ++i) {
    h[i] = (borders[i].second - borders[i].first) / static_cast<double>(n);
    cell_volume *= h[i];
  }

  int64_t total_points = 1;
  for (int i = 0; i < dim; ++i) {
    total_points *= static_cast<int64_t>(n);
  }

  int rank = 0;
  int world_size = 1;
  const bool is_mpi = ppc::util::IsUnderMpirun();

  if (is_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  }

  const int64_t chunk = total_points / world_size;
  const int64_t rem = total_points % world_size;
  const int64_t start = (rank * chunk) + std::min(rank, static_cast<int>(rem));
  const int64_t end = start + chunk + (rank < rem ? 1 : 0);

  double local_sum = 0.0;
  bool error_flag = false;

#pragma omp parallel for schedule(static) reduction(+ : local_sum) shared(error_flag) default(none) \
    firstprivate(func, borders, h, n, dim, start, end)
  for (int64_t linear_idx = start; linear_idx < end; ++linear_idx) {
    if (error_flag) {
      continue;
    }

    std::vector<double> point(dim);
    FillPoint(static_cast<size_t>(linear_idx), n, dim, borders, h, point);

    double f_val = func(point);
    if (!std::isfinite(f_val)) {
#pragma omp atomic write
      error_flag = true;
      continue;
    }
    local_sum += f_val;
  }

  if (error_flag) {
    return false;
  }

  double global_sum = local_sum;
  if (is_mpi && world_size > 1) {
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }

  GetOutput() = global_sum * cell_volume;
  return std::isfinite(GetOutput());
}

bool DergynovSIntegralsMultistepRectangleALL::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace dergynov_s_integrals_multistep_rectangle
