#include "shkrebko_m_calc_of_integral_rect/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
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

void ShkrebkoMCalcOfIntegralRectALL::DistributeData(int rank, size_t &dims) {
  if (!ppc::util::IsUnderMpirun()) {
    return;
  }

  int dims_int = static_cast<int>(dims);
  MPI_Bcast(&dims_int, 1, MPI_INT, 0, MPI_COMM_WORLD);
  dims = static_cast<size_t>(dims_int);

  if (rank != 0) {
    local_input_.limits.resize(dims);
    local_input_.n_steps.resize(dims);
  }

  for (size_t i = 0; i < dims; ++i) {
    MPI_Bcast(&local_input_.limits[i].first, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&local_input_.limits[i].second, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&local_input_.n_steps[i], 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
}

double ShkrebkoMCalcOfIntegralRectALL::ComputeChunkSum(size_t start_idx, size_t end_idx, const std::vector<double> &h,
                                                       const std::vector<std::pair<double, double>> &limits,
                                                       const std::vector<int> &n_steps,
                                                       const std::function<double(const std::vector<double> &)> &func) {
  if (start_idx >= end_idx) {
    return 0.0;
  }

  const size_t dim = limits.size();
  const size_t count = end_idx - start_idx;

  std::vector<int> indices(dim);
  std::vector<double> point(dim);

  size_t tmp = start_idx;
  for (int d = static_cast<int>(dim) - 1; d >= 0; --d) {
    indices[d] = static_cast<int>(tmp % static_cast<size_t>(n_steps[d]));
    tmp /= static_cast<size_t>(n_steps[d]);
  }

  double sum = 0.0;
  for (size_t iter = 0; iter < count; ++iter) {
    for (size_t d = 0; d < dim; ++d) {
      point[d] = limits[d].first + (static_cast<double>(indices[d]) + 0.5) * h[d];
    }
    sum += func(point);

    int d = static_cast<int>(dim) - 1;
    while (d >= 0) {
      indices[d]++;
      if (indices[d] < n_steps[d]) {
        break;
      }
      indices[d] = 0;
      d--;
    }
  }
  return sum;
}

bool ShkrebkoMCalcOfIntegralRectALL::RunImpl() {
  int rank = 0;
  int size = 1;
  bool is_mpi = ppc::util::IsUnderMpirun();

  if (is_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }

  size_t dims = (rank == 0) ? local_input_.limits.size() : 0;
  DistributeData(rank, dims);

  // защита от пустой функции
  if (!local_input_.func) {
    return false;
  }

  std::vector<double> h(dims);
  double cell_volume = 1.0;
  size_t total_points = 1;

  for (size_t i = 0; i < dims; ++i) {
    const double left = local_input_.limits[i].first;
    const double right = local_input_.limits[i].second;
    const int steps = local_input_.n_steps[i];
    h[i] = (right - left) / static_cast<double>(steps);
    cell_volume *= h[i];
    total_points *= static_cast<size_t>(steps);
  }

  size_t proc_chunk = total_points / static_cast<size_t>(size);
  size_t proc_remainder = total_points % static_cast<size_t>(size);
  size_t my_start = static_cast<size_t>(rank) * proc_chunk + std::min(static_cast<size_t>(rank), proc_remainder);
  size_t my_count = proc_chunk + (static_cast<size_t>(rank) < proc_remainder ? 1 : 0);

  double local_sum = 0.0;

  if (my_count > 0) {
    int num_threads = ppc::util::GetNumThreads();
    omp_set_num_threads(num_threads);

    const auto &func = local_input_.func;
    const auto &n_steps = local_input_.n_steps;
    const auto &limits = local_input_.limits;

#pragma omp parallel default(none) shared(h, dims, my_start, my_count, func, n_steps, limits) reduction(+ : local_sum)
    {
      int tid = omp_get_thread_num();
      int t_count = omp_get_num_threads();

      size_t thread_chunk = my_count / static_cast<size_t>(t_count);
      size_t thread_remainder = my_count % static_cast<size_t>(t_count);
      size_t thread_start =
          my_start + static_cast<size_t>(tid) * thread_chunk + std::min(static_cast<size_t>(tid), thread_remainder);
      size_t thread_count = thread_chunk + (static_cast<size_t>(tid) < thread_remainder ? 1 : 0);

      if (thread_count > 0) {
        local_sum += ComputeChunkSum(thread_start, thread_start + thread_count, h, limits, n_steps, func);
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
