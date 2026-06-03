#include "sannikov_i_integrals_rectangle_method/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>

#include "sannikov_i_integrals_rectangle_method/common/include/common.hpp"
#include "util/include/util.hpp"

namespace sannikov_i_integrals_rectangle_method {

SannikovIIntegralsRectangleMethodALL::SannikovIIntegralsRectangleMethodALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SannikovIIntegralsRectangleMethodALL::ValidationImpl() {
  const auto &[func, borders, num] = GetInput();
  if (borders.empty()) {
    return false;
  }
  for (const auto &[left_border, right_border] : borders) {
    if (!std::isfinite(left_border) || !std::isfinite(right_border)) {
      return false;
    }
    if (left_border >= right_border) {
      return false;
    }
  }

  return func && (num > 0) && (GetOutput() == 0.0);
}

bool SannikovIIntegralsRectangleMethodALL::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

namespace {

bool ComputeSlice(const std::function<double(const std::vector<double> &)> &funcx,
                  const std::vector<std::pair<double, double>> &brd, const std::vector<double> &h_vec,
                  std::size_t outer_dim, int64_t inner_cells, int num_splits, int outer, std::vector<int> &idx,
                  std::vector<double> &x, double &local_sum) {
  for (std::size_t i = 0; i < outer_dim; ++i) {
    idx[i] = 0;
  }
  idx[outer_dim] = outer;
  x[outer_dim] = brd[outer_dim].first + ((static_cast<double>(outer) + 0.5) * h_vec[outer_dim]);

  for (int64_t cell = 0; cell < inner_cells; ++cell) {
    for (std::size_t i = 0; i < outer_dim; ++i) {
      x[i] = brd[i].first + ((static_cast<double>(idx[i]) + 0.5) * h_vec[i]);
    }

    const double fx = funcx(x);
    if (!std::isfinite(fx)) {
      return false;
    }

    local_sum += fx;

    for (std::size_t pos = 0; pos < outer_dim; ++pos) {
      if (++idx[pos] < num_splits) {
        break;
      }
      idx[pos] = 0;
    }
  }
  return true;
}

}  // namespace

bool SannikovIIntegralsRectangleMethodALL::RunImpl() {
  const auto &[func, borders, num] = GetInput();
  const std::size_t dim = borders.size();

  const auto &funcx = func;
  const auto &brd = borders;
  const int num_splits = num;

  std::vector<double> h_vec(dim);
  double cell_v = 1.0;

  for (std::size_t i = 0; i < dim; ++i) {
    h_vec[i] = (brd[i].second - brd[i].first) / static_cast<double>(num_splits);
    if (!(h_vec[i] > 0.0) || !std::isfinite(h_vec[i])) {
      return false;
    }
    cell_v *= h_vec[i];
  }

  const std::size_t outer_dim = dim - 1;

  int64_t inner_cells = 1;
  for (std::size_t i = 0; i < outer_dim; ++i) {
    inner_cells *= num_splits;
  }

  int rank = 0;
  int num_procs = 1;
  const bool is_mpi = ppc::util::IsUnderMpirun();

  if (is_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  }

  const int chunk = num_splits / num_procs;
  const int rem = num_splits % num_procs;
  const int start = (rank * chunk) + std::min(rank, rem);
  const int end = start + chunk + (rank < rem ? 1 : 0);
  const int local_n = end - start;

  double local_sum = 0.0;
  bool error_flag = false;

  const int num_threads = ppc::util::GetNumThreads();
  const int max_threads = std::min(num_threads, local_n);
  std::vector<std::vector<int>> all_idx(max_threads, std::vector<int>(dim, 0));
  std::vector<std::vector<double>> all_x(max_threads, std::vector<double>(dim));

#pragma omp parallel for schedule(static) num_threads(max_threads) reduction(+ : local_sum) \
    shared(error_flag) default(none)                                                        \
    shared(funcx, brd, h_vec, outer_dim, inner_cells, num_splits, start, end, all_idx, all_x)
  for (int outer = start; outer < end; ++outer) {
    if (error_flag) {
      continue;
    }

    const int tid = omp_get_thread_num();
    double slice_sum = 0.0;

    if (!ComputeSlice(funcx, brd, h_vec, outer_dim, inner_cells, num_splits, outer, all_idx[tid], all_x[tid],
                      slice_sum)) {
      error_flag = true;
    }

    local_sum += slice_sum;
  }

  double global_sum = local_sum;
  if (is_mpi) {
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }

  if (error_flag) {
    return false;
  }

  GetOutput() = global_sum * cell_v;
  return std::isfinite(GetOutput());
}

bool SannikovIIntegralsRectangleMethodALL::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace sannikov_i_integrals_rectangle_method
