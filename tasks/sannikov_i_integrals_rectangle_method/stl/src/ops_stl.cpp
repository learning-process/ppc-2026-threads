#include "sannikov_i_integrals_rectangle_method/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <thread>
#include <utility>
#include <vector>

#include "sannikov_i_integrals_rectangle_method/common/include/common.hpp"

namespace sannikov_i_integrals_rectangle_method {

SannikovIIntegralsRectangleMethodSTL::SannikovIIntegralsRectangleMethodSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SannikovIIntegralsRectangleMethodSTL::ValidationImpl() {
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

bool SannikovIIntegralsRectangleMethodSTL::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

namespace {

bool ComputeSlice(const std::function<double(const std::vector<double> &)> &funcx,
                  const std::vector<std::pair<double, double>> &brd, const std::vector<double> &sz,
                  std::size_t outer_dim, int64_t inner_cells, int num_splits, int outer, std::vector<int> &idx,
                  std::vector<double> &x_vec, double &local_sum) {
  for (std::size_t i = 0; i < outer_dim; ++i) {
    idx[i] = 0;
  }
  idx[outer_dim] = outer;
  x_vec[outer_dim] = brd[outer_dim].first + ((static_cast<double>(outer) + 0.5) * sz[outer_dim]);

  for (int64_t cell = 0; cell < inner_cells; ++cell) {
    for (std::size_t i = 0; i < outer_dim; ++i) {
      x_vec[i] = brd[i].first + ((static_cast<double>(idx[i]) + 0.5) * sz[i]);
    }

    const double func_x = funcx(x_vec);
    if (!std::isfinite(func_x)) {
      return false;
    }

    local_sum += func_x;

    for (std::size_t position = 0; position < outer_dim; ++position) {
      if (++idx[position] < num_splits) {
        break;
      }
      idx[position] = 0;
    }
  }
  return true;
}

void ThreadWorker(const std::function<double(const std::vector<double> &)> &funcx,
                  const std::vector<std::pair<double, double>> &brd, const std::vector<double> &sz, std::size_t dim,
                  std::size_t outer_dim, int64_t inner_cells, int num_splits, int start, int end, double &partial_sum,
                  int &error_flag) {
  std::vector<int> idx(dim, 0);
  std::vector<double> x_vec(dim);
  double thread_sum = 0.0;

  for (int outer = start; outer < end; ++outer) {
    double local_sum = 0.0;
    if (!ComputeSlice(funcx, brd, sz, outer_dim, inner_cells, num_splits, outer, idx, x_vec, local_sum)) {
      error_flag = 1;
      return;
    }
    thread_sum += local_sum;
  }

  partial_sum = thread_sum;
}

}  // namespace

bool SannikovIIntegralsRectangleMethodSTL::RunImpl() {
  const auto &[func, borders, num] = GetInput();
  const std::size_t dim = borders.size();

  const auto &funcx = func;
  const auto &brd = borders;
  const int num_splits = num;

  std::vector<double> sz(dim);
  double cell_v = 1.0;

  for (std::size_t i = 0; i < dim; ++i) {
    sz[i] = (brd[i].second - brd[i].first) / static_cast<double>(num_splits);
    if (!(sz[i] > 0.0) || !std::isfinite(sz[i])) {
      return false;
    }
    cell_v *= sz[i];
  }

  const std::size_t outer_dim = dim - 1;

  int64_t inner_cells = 1;
  for (std::size_t i = 0; i < outer_dim; ++i) {
    inner_cells *= num_splits;
  }

  const int num_threads = static_cast<int>(std::thread::hardware_concurrency());
  const int actual_threads = std::min(num_threads > 0 ? num_threads : 4, num_splits);

  std::vector<double> partial_sums(actual_threads, 0.0);
  std::vector<int> error_flags(actual_threads, 0);
  std::vector<std::thread> threads(actual_threads);

  for (int tid = 0; tid < actual_threads; ++tid) {
    const int chunk = num_splits / actual_threads;
    const int rem = num_splits % actual_threads;
    const int start = (tid * chunk) + std::min(tid, rem);
    const int end = start + chunk + (tid < rem ? 1 : 0);

    threads[tid] = std::thread([&funcx, &brd, &sz, dim, outer_dim, inner_cells, num_splits, start, end,
                                &partial_sum = partial_sums[tid], &error_flag = error_flags[tid]]() {
      ThreadWorker(funcx, brd, sz, dim, outer_dim, inner_cells, num_splits, start, end, partial_sum, error_flag);
    });
  }

  for (auto &thr : threads) {
    thr.join();
  }

  for (int tid = 0; tid < actual_threads; ++tid) {
    if (error_flags[tid] != 0) {
      return false;
    }
  }

  double sum = 0.0;
  for (int tid = 0; tid < actual_threads; ++tid) {
    sum += partial_sums[tid];
  }

  GetOutput() = sum * cell_v;
  return std::isfinite(GetOutput());
}

bool SannikovIIntegralsRectangleMethodSTL::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace sannikov_i_integrals_rectangle_method
