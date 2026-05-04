#include "sannikov_i_integrals_rectangle_method/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>

#include "sannikov_i_integrals_rectangle_method/common/include/common.hpp"

namespace sannikov_i_integrals_rectangle_method {

SannikovIIntegralsRectangleMethodOMP::SannikovIIntegralsRectangleMethodOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SannikovIIntegralsRectangleMethodOMP::ValidationImpl() {
  const auto &[func, borders, n] = GetInput();
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

  return func && (n > 0) && (GetOutput() == 0.0);
}

bool SannikovIIntegralsRectangleMethodOMP::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

namespace {

bool ComputeSlice(const std::function<double(const std::vector<double> &)> &f,
                  const std::vector<std::pair<double, double>> &brd, const std::vector<double> &h,
                  std::size_t outer_dim, int64_t inner_cells, int num_splits, int outer, std::vector<int> &idx,
                  std::vector<double> &x, double &local_sum) {
  for (std::size_t i = 0; i < outer_dim; ++i) {
    idx[i] = 0;
  }
  idx[outer_dim] = outer;
  x[outer_dim] = brd[outer_dim].first + ((static_cast<double>(outer) + 0.5) * h[outer_dim]);

  for (int64_t cell = 0; cell < inner_cells; ++cell) {
    for (std::size_t i = 0; i < outer_dim; ++i) {
      x[i] = brd[i].first + ((static_cast<double>(idx[i]) + 0.5) * h[i]);
    }

    const double fx = f(x);
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

bool SannikovIIntegralsRectangleMethodOMP::RunImpl() {
  const auto &[func, borders, n] = GetInput();
  const std::size_t dim = borders.size();

  const auto &f = func;
  const auto &brd = borders;
  const int num_splits = n;

  std::vector<double> h(dim);
  double cell_v = 1.0;

  for (std::size_t i = 0; i < dim; ++i) {
    h[i] = (brd[i].second - brd[i].first) / static_cast<double>(num_splits);
    if (!(h[i] > 0.0) || !std::isfinite(h[i])) {
      return false;
    }
    cell_v *= h[i];
  }

  const std::size_t outer_dim = dim - 1;

  int64_t inner_cells = 1;
  for (std::size_t i = 0; i < outer_dim; ++i) {
    inner_cells *= num_splits;
  }

  const int max_threads = omp_get_max_threads();
  std::vector<std::vector<int>> all_idx(max_threads, std::vector<int>(dim, 0));
  std::vector<std::vector<double>> all_x(max_threads, std::vector<double>(dim));

  double sum = 0.0;
  bool error_flag = false;

#pragma omp parallel for schedule(static) reduction(+ : sum) shared(error_flag) default(none) \
    shared(f, brd, h, outer_dim, inner_cells, num_splits, all_idx, all_x)
  for (int outer = 0; outer < num_splits; ++outer) {
    if (error_flag) {
      continue;
    }

    const int tid = omp_get_thread_num();
    double local_sum = 0.0;

    if (!ComputeSlice(f, brd, h, outer_dim, inner_cells, num_splits, outer, all_idx[tid], all_x[tid], local_sum)) {
      error_flag = true;
    }

    sum += local_sum;
  }

  if (error_flag) {
    return false;
  }

  GetOutput() = sum * cell_v;
  return std::isfinite(GetOutput());
}

bool SannikovIIntegralsRectangleMethodOMP::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace sannikov_i_integrals_rectangle_method
