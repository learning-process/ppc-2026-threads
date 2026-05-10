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

  use_mpi_ = ppc::util::IsUnderMpirun();
  if (use_mpi_) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size_);
  } else {
    rank_ = 0;
    world_size_ = 1;
  }

  return true;
}

void ShkrebkoMCalcOfIntegralRectALL::BroadcastInputData() {
  int num_dims = (rank_ == 0) ? static_cast<int>(local_input_.limits.size()) : 0;
  MPI_Bcast(&num_dims, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank_ != 0) {
    local_input_.limits.resize(static_cast<std::size_t>(num_dims));
    local_input_.n_steps.resize(static_cast<std::size_t>(num_dims));
  }

  const auto udims = static_cast<std::size_t>(num_dims);
  std::vector<double> flat_limits(udims * 2);

  if (rank_ == 0) {
    for (std::size_t i = 0; i < udims; ++i) {
      flat_limits[(2 * i)] = local_input_.limits[i].first;
      flat_limits[(2 * i) + 1] = local_input_.limits[i].second;
    }
  }

  MPI_Bcast(flat_limits.data(), 2 * num_dims, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(local_input_.n_steps.data(), num_dims, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank_ != 0) {
    for (std::size_t i = 0; i < udims; ++i) {
      local_input_.limits[i].first = flat_limits[(2 * i)];
      local_input_.limits[i].second = flat_limits[(2 * i) + 1];
    }
  }
}

bool ShkrebkoMCalcOfIntegralRectALL::RunImpl() {
  if (use_mpi_) {
    auto saved_func = local_input_.func;
    BroadcastInputData();
    if (!local_input_.func) {
      local_input_.func = saved_func;
    }
  }

  const std::size_t dims = local_input_.limits.size();
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

  const std::int64_t base_per_rank = total_points / world_size_;
  const std::int64_t leftover = total_points % world_size_;

  const std::int64_t my_start = (rank_ * base_per_rank) + std::min(static_cast<std::int64_t>(rank_), leftover);
  const std::int64_t my_count = base_per_rank + (rank_ < static_cast<int>(leftover) ? 1 : 0);

  double rank_sum = 0.0;

  if (my_count > 0) {
    omp_set_num_threads(ppc::util::GetNumThreads());

#pragma omp parallel default(none) shared(my_start, my_count, h, limits, n_steps, func, dims) reduction(+ : rank_sum)
    {
      std::vector<double> point(dims);

#pragma omp for schedule(static)
      for (std::int64_t idx = 0; idx < my_count; ++idx) {
        auto global_idx = static_cast<std::size_t>(my_start + idx);

        std::size_t tmp = global_idx;
        for (int dim = static_cast<int>(dims) - 1; dim >= 0; --dim) {
          auto sc = static_cast<std::size_t>(n_steps[dim]);
          const std::size_t ci = tmp % sc;
          tmp /= sc;
          point[dim] = limits[dim].first + ((static_cast<double>(ci) + 0.5) * h[dim]);
        }

        rank_sum += func(point);
      }
    }
  }

  if (use_mpi_) {
    double aggregated_sum = 0.0;
    MPI_Reduce(&rank_sum, &aggregated_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank_ == 0) {
      res_ = aggregated_sum * cell_volume;
    }
  } else {
    res_ = rank_sum * cell_volume;
  }

  return true;
}

bool ShkrebkoMCalcOfIntegralRectALL::PostProcessingImpl() {
  if (rank_ == 0) {
    GetOutput() = res_;
  }
  return true;
}

}  // namespace shkrebko_m_calc_of_integral_rect
