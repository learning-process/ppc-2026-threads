#include "shkrebko_m_calc_of_integral_rect/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
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
  int rank = 0;
  if (ppc::util::IsUnderMpirun()) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }

  if (rank != 0) {
    return true;
  }

  const auto &input = GetInput();

  if (!input.func) {
    return false;
  }
  if (input.limits.size() != input.n_steps.size() || input.limits.empty()) {
    return false;
  }
  if (!std::ranges::all_of(input.n_steps, [](int n) { return n > 0; })) {
    return false;
  }
  if (!std::ranges::all_of(input.limits, [](const auto &lim) { return lim.first < lim.second; })) {
    return false;
  }

  return true;
}

bool ShkrebkoMCalcOfIntegralRectALL::PreProcessingImpl() {
  local_input_ = GetInput();
  res_ = 0.0;
  return true;
}

void ShkrebkoMCalcOfIntegralRectALL::BroadcastInputData(int rank, std::size_t &dims) {
  if (!ppc::util::IsUnderMpirun()) {
    return;
  }

  int dims_io = static_cast<int>(dims);
  MPI_Bcast(&dims_io, 1, MPI_INT, 0, MPI_COMM_WORLD);
  dims = static_cast<std::size_t>(dims_io);

  if (rank != 0) {
    local_input_.limits.resize(dims);
    local_input_.n_steps.resize(dims);
  }

  std::vector<double> flat_limits(dims * 2);
  if (rank == 0) {
    for (std::size_t i = 0; i < dims; ++i) {
      flat_limits[2 * i] = local_input_.limits[i].first;
      flat_limits[(2 * i) + 1] = local_input_.limits[i].second;
    }
  }

  MPI_Bcast(flat_limits.data(), dims_io * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(local_input_.n_steps.data(), dims_io, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    for (std::size_t i = 0; i < dims; ++i) {
      local_input_.limits[i].first = flat_limits[2 * i];
      local_input_.limits[i].second = flat_limits[(2 * i) + 1];
    }
  }
}

void ShkrebkoMCalcOfIntegralRectALL::FlatIndexToPoint(std::size_t flat_idx, const std::vector<double> &h,
                                                      std::vector<double> &point) const {
  const std::size_t dims = local_input_.limits.size();
  std::size_t remainder = flat_idx;

  for (int d = static_cast<int>(dims) - 1; d >= 0; --d) {
    const std::size_t step_count = static_cast<std::size_t>(local_input_.n_steps[d]);
    const std::size_t coord = remainder % step_count;
    remainder /= step_count;
    point[d] = local_input_.limits[d].first + ((static_cast<double>(coord) + 0.5) * h[d]);
  }
}

bool ShkrebkoMCalcOfIntegralRectALL::RunImpl() {
  int rank = 0;
  int size = 1;
  const bool is_mpi = ppc::util::IsUnderMpirun();

  if (is_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }

  std::size_t dims = (rank == 0) ? local_input_.limits.size() : 0;
  BroadcastInputData(rank, dims);

  if (!local_input_.func) {
    return false;
  }

  const auto &limits = local_input_.limits;
  const auto &n_steps = local_input_.n_steps;
  const auto &func = local_input_.func;

  std::vector<double> h(dims);
  double cell_volume = 1.0;
  std::int64_t total_points = 1;

  for (std::size_t i = 0; i < dims; ++i) {
    h[i] = (limits[i].second - limits[i].first) / static_cast<double>(n_steps[i]);
    cell_volume *= h[i];
    total_points *= static_cast<std::int64_t>(n_steps[i]);
  }

  const std::int64_t base_per_rank = total_points / static_cast<std::int64_t>(size);
  const std::int64_t leftover = total_points % static_cast<std::int64_t>(size);

  const std::int64_t my_start =
      (static_cast<std::int64_t>(rank) * base_per_rank) + std::min<std::int64_t>(rank, leftover);
  const std::int64_t my_count = base_per_rank + ((static_cast<std::int64_t>(rank) < leftover) ? 1 : 0);

  double local_sum = 0.0;

  if (my_count > 0) {
    omp_set_num_threads(ppc::util::GetNumThreads());

#pragma omp parallel default(none) shared(my_start, my_count, h, dims, func) reduction(+ : local_sum)
    {
      std::vector<double> point(dims);

#pragma omp for schedule(static)
      for (std::int64_t idx = 0; idx < my_count; ++idx) {
        const std::size_t global_idx = static_cast<std::size_t>(my_start + idx);
        FlatIndexToPoint(global_idx, h, point);
        local_sum += func(point);
      }
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
