#include "../include/rect_method_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

#include "../../common/include/common.hpp"
#include "util/include/util.hpp"

namespace kutergin_v_multidimensional_integration_rect_method {

RectMethodALL::RectMethodALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool RectMethodALL::ValidationImpl() {
  const auto &input = GetInput();
  if (input.limits.size() != input.n_steps.size() || input.limits.empty()) {
    return false;
  }
  return std::ranges::all_of(input.n_steps, [](int n) { return n > 0; });
}

bool RectMethodALL::PreProcessingImpl() {
  local_input_ = GetInput();
  res_ = 0.0;
  return true;
}

bool RectMethodALL::RunImpl() {
  int rank = 0;
  int size = 1;
  bool is_mpi = ppc::util::IsUnderMpirun();

  if (is_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }

  size_t dims = (rank == 0) ? local_input_.limits.size() : 0;

  DistributeData(rank, dims);

  size_t total_iterations = 1;
  std::vector<double> h(dims);
  double d_v = 1.0;

  for (size_t i = 0; i < dims; ++i) {
    total_iterations *= local_input_.n_steps[i];
    h[i] = (local_input_.limits[i].second - local_input_.limits[i].first) / local_input_.n_steps[i];
    d_v *= h[i];
  }

  size_t proc_chunk = total_iterations / size;
  size_t proc_remainder = total_iterations % size;

  size_t my_proc_count = proc_chunk + (std::cmp_less(rank, proc_remainder) ? 1 : 0);
  size_t my_proc_start = (rank * proc_chunk) + std::min(static_cast<size_t>(rank), proc_remainder);

  double local_sum = 0.0;

  if (my_proc_count > 0) {
    int num_threads = ppc::util::GetNumThreads();
    omp_set_num_threads(num_threads);

    const auto &func = local_input_.func;
    const auto &n_steps = local_input_.n_steps;
    const auto &limits = local_input_.limits;

#pragma omp parallel default(none) shared(h, dims, my_proc_start, my_proc_count, func, n_steps, limits) \
    reduction(+ : local_sum)
    {
      int tid = omp_get_thread_num();
      int t_count = omp_get_num_threads();

      size_t thread_chunk = my_proc_count / t_count;
      size_t thread_remainder = my_proc_count % t_count;

      size_t my_thread_count = thread_chunk + (std::cmp_less(tid, thread_remainder) ? 1 : 0);
      size_t my_thread_start =
          my_proc_start + (tid * thread_chunk) + std::min(static_cast<size_t>(tid), thread_remainder);

      local_sum += CalculateChunkSum(my_thread_start, my_thread_start + my_thread_count, h, limits, n_steps, func);
    }
  }

  if (is_mpi) {
    double global_sum = 0.0;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
      res_ = global_sum * d_v;
    }
  } else {
    res_ = local_sum * d_v;
  }

  return true;
}

bool RectMethodALL::PostProcessingImpl() {
  int rank = 0;
  if (ppc::util::IsUnderMpirun()) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }
  if (rank == 0) {
    GetOutput() = res_;
  }
  return true;
}

void RectMethodALL::DistributeData(int rank, size_t &dims) {
  if (!ppc::util::IsUnderMpirun()) {
    return;
  }

  int dims_io = static_cast<int>(dims);
  MPI_Bcast(&dims_io, 1, MPI_INT, 0, MPI_COMM_WORLD);
  dims = static_cast<size_t>(dims_io);

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

double RectMethodALL::CalculateChunkSum(size_t start_idx, size_t end_idx, const std::vector<double> &h,
                                        const std::vector<std::pair<double, double>> &limits,
                                        const std::vector<int> &n_steps,
                                        const std::function<double(const std::vector<double> &)> &func) {
  if (start_idx >= end_idx) {
    return 0.0;
  }

  size_t count = end_idx - start_idx;
  size_t dims = limits.size();  // число размерностей пространства
  std::vector<int> current_indices(dims, 0);
  std::vector<double> coords(dims);  // создание вектора координат размером dims

  size_t temp_idx = start_idx;
  for (int dim = static_cast<int>(dims) - 1; dim >= 0; --dim) {
    current_indices[dim] = static_cast<int>(temp_idx % n_steps[dim]);
    temp_idx /= n_steps[dim];
  }

  double chunk_sum = 0.0;

  for (size_t i = 0; i < count; ++i) {
    for (size_t dm = 0; dm < dims; ++dm) {
      coords[dm] = limits[dm].first + ((current_indices[dm] + 0.5) * h[dm]);  // реальная координата
    }

    chunk_sum += func(coords);  // вычисление функции в точке

    for (int dm = static_cast<int>(dims) - 1; dm >= 0; --dm) {
      current_indices[dm]++;
      if (current_indices[dm] < n_steps[dm]) {
        break;
      }
      current_indices[dm] = 0;
    }
  }

  return chunk_sum;
}

}  // namespace kutergin_v_multidimensional_integration_rect_method
