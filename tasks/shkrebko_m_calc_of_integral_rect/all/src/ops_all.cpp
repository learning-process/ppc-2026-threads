#include "shkrebko_m_calc_of_integral_rect/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
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

  std::vector<double> left_bounds(dims);
  std::vector<double> right_bounds(dims);

  if (rank == 0) {
    for (std::size_t i = 0; i < dims; ++i) {
      left_bounds[i] = local_input_.limits[i].first;
      right_bounds[i] = local_input_.limits[i].second;
    }
  }

  MPI_Bcast(left_bounds.data(), dims_io, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(right_bounds.data(), dims_io, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(local_input_.n_steps.data(), dims_io, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    for (std::size_t i = 0; i < dims; ++i) {
      local_input_.limits[i].first = left_bounds[i];
      local_input_.limits[i].second = right_bounds[i];
    }
  }
}

std::size_t ShkrebkoMCalcOfIntegralRectALL::SelectSplitDimension() const {
  const auto &n_steps = local_input_.n_steps;
  const auto it = std::ranges::max_element(n_steps);
  return static_cast<std::size_t>(std::distance(n_steps.begin(), it));
}

bool ShkrebkoMCalcOfIntegralRectALL::ComputeSliceSum(std::size_t fixed_dim, std::size_t fixed_idx,
                                                     const std::vector<double> &h, double &slice_sum) const {
  slice_sum = 0.0;

  const std::size_t dim = local_input_.limits.size();
  const auto &limits = local_input_.limits;
  const auto &n_steps = local_input_.n_steps;
  const auto &func = local_input_.func;

  if (!func) {
    return false;
  }

  std::vector<double> point(dim);
  point[fixed_dim] = limits[fixed_dim].first + ((static_cast<double>(fixed_idx) + 0.5) * h[fixed_dim]);

  if (dim == 1) {
    double f_val = func(point);
    if (!std::isfinite(f_val)) {
      return false;
    }
    slice_sum = f_val;
    return true;
  }

  std::vector<std::size_t> free_dims;
  free_dims.reserve(dim - 1);
  for (std::size_t i = 0; i < dim; ++i) {
    if (i != fixed_dim) {
      free_dims.push_back(i);
    }
  }

  const auto free_count = free_dims.size();
  std::vector<int> idx(free_dims.size(), 0);

  while (true) {
    for (std::size_t j = 0; j < free_dims.size(); ++j) {
      std::size_t d = free_dims[j];
      point[d] = limits[d].first + ((static_cast<double>(idx[j]) + 0.5) * h[d]);
    }

    double f_val = func(point);
    if (!std::isfinite(f_val)) {
      return false;
    }
    slice_sum += f_val;

    // Одометр – инкремент комбинации
    int level = static_cast<int>(free_count) - 1;
    while (level >= 0) {
      if (++idx[level] < n_steps[free_dims[static_cast<std::size_t>(level)]]) {
        break;
      }
      idx[level] = 0;
      --level;
    }
    if (level < 0) {
      break;
    }
  }

  return true;
}

double ShkrebkoMCalcOfIntegralRectALL::ComputeCellVolume(const std::vector<double> &h) const {
  double volume = 1.0;
  for (double step : h) {
    volume *= step;
  }
  return volume;
}

std::vector<std::size_t> ShkrebkoMCalcOfIntegralRectALL::DistributeSlices(int rank, int size,
                                                                          std::size_t split_steps) const {
  std::vector<std::size_t> local_slices;
  for (std::size_t slice = static_cast<std::size_t>(rank); slice < split_steps;
       slice += static_cast<std::size_t>(size)) {
    local_slices.push_back(slice);
  }
  return local_slices;
}

std::pair<double, bool> ShkrebkoMCalcOfIntegralRectALL::ComputeLocalSum(const std::vector<std::size_t> &local_slices,
                                                                        const std::vector<double> &h,
                                                                        std::size_t split_dim) const {
  double local_sum = 0.0;
  bool local_ok = true;

  if (local_slices.empty()) {
    return {local_sum, local_ok};
  }

  const auto requested = static_cast<std::size_t>(std::max(1, ppc::util::GetNumThreads()));
  const std::size_t thread_count = std::min(requested, local_slices.size());

  std::vector<std::thread> threads;
  threads.reserve(thread_count);
  std::vector<double> partial_sums(thread_count, 0.0);
  std::vector<bool> partial_ok(thread_count, true);

  for (std::size_t tid = 0; tid < thread_count; ++tid) {
    threads.emplace_back([this, &local_slices, &h, &partial_sums, &partial_ok, tid, thread_count, split_dim]() {
      double thread_sum = 0.0;
      for (std::size_t pos = tid; pos < local_slices.size(); pos += thread_count) {
        double slice_sum = 0.0;
        if (!ComputeSliceSum(split_dim, local_slices[pos], h, slice_sum)) {
          partial_ok[tid] = false;
          return;
        }
        thread_sum += slice_sum;
      }
      partial_sums[tid] = thread_sum;
    });
  }

  for (auto &t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }

  for (std::size_t i = 0; i < thread_count; ++i) {
    local_sum += partial_sums[i];
    if (!partial_ok[i]) {
      local_ok = false;
    }
  }

  return {local_sum, local_ok};
}

bool ShkrebkoMCalcOfIntegralRectALL::FinalizeResult(double local_sum, double cell_volume, bool local_ok, bool is_mpi,
                                                    int rank) {
  if (is_mpi) {
    int local_ok_int = local_ok ? 1 : 0;
    int global_ok = 0;
    MPI_Allreduce(&local_ok_int, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (global_ok == 0) {
      return false;
    }

    double global_sum = 0.0;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
      res_ = global_sum * cell_volume;
    }
  } else {
    if (!local_ok) {
      return false;
    }
    res_ = local_sum * cell_volume;
  }
  return true;
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

  const auto &limits = local_input_.limits;
  const auto &n_steps = local_input_.n_steps;

  std::vector<double> h(dims);
  for (std::size_t i = 0; i < dims; ++i) {
    h[i] = (limits[i].second - limits[i].first) / static_cast<double>(n_steps[i]);
  }
  const double cell_volume = ComputeCellVolume(h);

  const std::size_t split_dim = SelectSplitDimension();
  const std::size_t split_steps = static_cast<std::size_t>(n_steps[split_dim]);
  const auto local_slices = DistributeSlices(rank, size, split_steps);

  const auto [local_sum, local_ok] = ComputeLocalSum(local_slices, h, split_dim);
  return FinalizeResult(local_sum, cell_volume, local_ok, is_mpi, rank);
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
