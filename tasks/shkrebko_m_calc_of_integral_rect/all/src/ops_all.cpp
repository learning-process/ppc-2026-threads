#include "shkrebko_m_calc_of_integral_rect/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <thread>
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
  if (!use_mpi_) {
    return;
  }

  int dims_int = (rank_ == 0) ? static_cast<int>(local_input_.limits.size()) : 0;
  MPI_Bcast(&dims_int, 1, MPI_INT, 0, MPI_COMM_WORLD);

  const std::size_t dims = static_cast<std::size_t>(dims_int);

  if (rank_ != 0) {
    local_input_.limits.resize(dims);
    local_input_.n_steps.resize(dims);
  }

  std::vector<double> left_bounds(dims);
  std::vector<double> right_bounds(dims);

  if (rank_ == 0) {
    for (std::size_t i = 0; i < dims; ++i) {
      left_bounds[i] = local_input_.limits[i].first;
      right_bounds[i] = local_input_.limits[i].second;
    }
  }

  MPI_Bcast(left_bounds.data(), dims_int, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(right_bounds.data(), dims_int, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(local_input_.n_steps.data(), dims_int, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank_ != 0) {
    for (std::size_t i = 0; i < dims; ++i) {
      local_input_.limits[i] = {left_bounds[i], right_bounds[i]};
    }
  }
}

std::size_t ShkrebkoMCalcOfIntegralRectALL::SelectSplitDimension() const {
  const auto &n_steps = local_input_.n_steps;
  const auto it = std::max_element(n_steps.begin(), n_steps.end());
  return static_cast<std::size_t>(std::distance(n_steps.begin(), it));
}

double ShkrebkoMCalcOfIntegralRectALL::ComputeSliceSum(std::size_t fixed_dim, std::size_t fixed_idx,
                                                       const std::vector<double> &h) const {
  const std::size_t dim = local_input_.limits.size();
  const auto &limits = local_input_.limits;
  const auto &n_steps = local_input_.n_steps;
  const auto &func = local_input_.func;

  std::vector<double> point(dim);
  point[fixed_dim] = limits[fixed_dim].first + ((static_cast<double>(fixed_idx) + 0.5) * h[fixed_dim]);

  if (dim == 1) {
    return func(point);
  }

  std::vector<int> idx(dim, 0);
  double sum = 0.0;

  while (true) {
    for (std::size_t i = 0; i < dim; ++i) {
      if (i == fixed_dim) {
        continue;
      }
      point[i] = limits[i].first + ((static_cast<double>(idx[i]) + 0.5) * h[i]);
    }

    sum += func(point);

    int level = static_cast<int>(dim) - 1;
    while (level >= 0) {
      if (static_cast<std::size_t>(level) == fixed_dim) {
        level--;
        continue;
      }

      idx[level]++;
      if (idx[level] < n_steps[level]) {
        break;
      }

      idx[level] = 0;
      level--;
    }

    if (level < 0) {
      break;
    }
  }

  return sum;
}

bool ShkrebkoMCalcOfIntegralRectALL::RunImpl() {
  if (use_mpi_) {
    BroadcastInputData();
  }

  const std::size_t dim = local_input_.limits.size();
  const auto &limits = local_input_.limits;
  const auto &n_steps = local_input_.n_steps;

  std::vector<double> h(dim);
  double cell_volume = 1.0;
  for (std::size_t i = 0; i < dim; ++i) {
    h[i] = (limits[i].second - limits[i].first) / static_cast<double>(n_steps[i]);
    cell_volume *= h[i];
  }

  const std::size_t split_dim = SelectSplitDimension();
  const std::size_t split_steps = static_cast<std::size_t>(n_steps[split_dim]);

  const std::size_t requested_threads = static_cast<std::size_t>(std::max(1, ppc::util::GetNumThreads()));
  const std::size_t thread_count = std::min(requested_threads, std::max<std::size_t>(1, split_steps));

  std::vector<std::thread> threads;
  threads.reserve(thread_count);

  std::vector<double> partial_sums(thread_count, 0.0);

  const std::size_t worker_stride = static_cast<std::size_t>(world_size_) * thread_count;

  for (std::size_t tid = 0; tid < thread_count; ++tid) {
    threads.emplace_back([&, tid, split_dim, worker_stride, split_steps]() {
      const std::size_t worker_offset = (static_cast<std::size_t>(rank_) * thread_count) + tid;

      double local_sum = 0.0;
      for (std::size_t fixed_idx = worker_offset; fixed_idx < split_steps; fixed_idx += worker_stride) {
        local_sum += ComputeSliceSum(split_dim, fixed_idx, h);
      }

      partial_sums[tid] = local_sum;
    });
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }

  double local_sum = 0.0;
  for (double val : partial_sums) {
    local_sum += val;
  }

  if (use_mpi_) {
    double global_sum = 0.0;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank_ == 0) {
      res_ = global_sum * cell_volume;
    }
  } else {
    res_ = local_sum * cell_volume;
  }

  return true;
}

bool ShkrebkoMCalcOfIntegralRectALL::PostProcessingImpl() {
  if (!use_mpi_ || rank_ == 0) {
    GetOutput() = res_;
  }
  return true;
}

}  // namespace shkrebko_m_calc_of_integral_rect
