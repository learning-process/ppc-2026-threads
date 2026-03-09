#include "galkin_d_multidim_integrals_rectangles/omp/include/ops_omp.hpp"

#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

#include "galkin_d_multidim_integrals_rectangles/common/include/common.hpp"
#include "util/include/util.hpp"

namespace galkin_d_multidim_integrals_rectangles {

GalkinDMultidimIntegralsRectanglesOMP::GalkinDMultidimIntegralsRectanglesOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool GalkinDMultidimIntegralsRectanglesOMP::ValidationImpl() {
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

bool GalkinDMultidimIntegralsRectanglesOMP::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool GalkinDMultidimIntegralsRectanglesOMP::RunImpl() {
  const InType &input = GetInput();
  const auto &func = std::get<0>(input);
  const auto &borders = std::get<1>(input);
  const int n = std::get<2>(input);
  const std::size_t dim = borders.size();

  std::vector<double> h(dim);
  double cell_v = 1.0;

  for (std::size_t i = 0; i < dim; ++i) {
    const double left_border = borders[i].first;
    const double right_border = borders[i].second;

    h[i] = (right_border - left_border) / static_cast<double>(n);
    if (!(h[i] > 0.0) || !std::isfinite(h[i])) {
      return false;
    }

    cell_v *= h[i];
  }

  std::size_t total_cells = 1;
  for (std::size_t i = 0; i < dim; ++i) {
    if (total_cells > (std::numeric_limits<std::size_t>::max() / static_cast<std::size_t>(n))) {
      return false;
    }
    total_cells *= static_cast<std::size_t>(n);
  }
  if (total_cells > static_cast<std::size_t>(LLONG_MAX)) {
    return false;
  }

  double sum = 0.0;
  const int requested_threads = ppc::util::GetNumThreads();
  const int positive_threads = (requested_threads > 0) ? requested_threads : 1;
  const int max_threads_by_work = (total_cells > static_cast<std::size_t>(std::numeric_limits<int>::max()))
                                      ? std::numeric_limits<int>::max()
                                      : static_cast<int>(total_cells);
  const int num_threads =
      std::cmp_greater(positive_threads, max_threads_by_work) ? max_threads_by_work : positive_threads;
  if (num_threads <= 0) {
    return false;
  }
  const auto total_cells_i64 = static_cast<std::int64_t>(total_cells);

#pragma omp parallel for default(none) shared(borders, h, dim, func, n, total_cells_i64, num_threads) \
    reduction(+ : sum) schedule(static) num_threads(num_threads)
  for (std::int64_t linear_idx = 0; linear_idx < total_cells_i64; ++linear_idx) {
    std::vector<double> x(dim);
    auto tmp = static_cast<std::size_t>(linear_idx);

    for (std::size_t i = 0; i < dim; ++i) {
      const std::size_t idx_i = tmp % static_cast<std::size_t>(n);
      tmp /= static_cast<std::size_t>(n);
      x[i] = borders[i].first + ((static_cast<double>(idx_i) + 0.5) * h[i]);
    }

    sum += func(x);
  }

  GetOutput() = sum * cell_v;

  return std::isfinite(GetOutput());
}

bool GalkinDMultidimIntegralsRectanglesOMP::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace galkin_d_multidim_integrals_rectangles
