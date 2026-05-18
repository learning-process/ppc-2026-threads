#include "galkin_d_multidim_integrals_rectangles/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

#include "galkin_d_multidim_integrals_rectangles/common/include/common.hpp"

namespace galkin_d_multidim_integrals_rectangles {

GalkinDMultidimIntegralsRectanglesTBB::GalkinDMultidimIntegralsRectanglesTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool GalkinDMultidimIntegralsRectanglesTBB::ValidationImpl() {
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

bool GalkinDMultidimIntegralsRectanglesTBB::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool GalkinDMultidimIntegralsRectanglesTBB::RunImpl() {
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

  const auto total_cells_i64 = static_cast<std::int64_t>(total_cells);

  double sum = tbb::parallel_reduce(tbb::blocked_range<std::int64_t>(0, total_cells_i64), 0.0,
                                    [&](const tbb::blocked_range<std::int64_t> &r, double local_sum) {
    std::vector<double> x(dim);

    for (std::int64_t linear_idx = r.begin(); linear_idx != r.end(); ++linear_idx) {
      auto tmp = static_cast<std::size_t>(linear_idx);

      for (std::size_t i = 0; i < dim; ++i) {
        const std::size_t idx_i = tmp % static_cast<std::size_t>(n);
        tmp /= static_cast<std::size_t>(n);
        x[i] = borders[i].first + ((static_cast<double>(idx_i) + 0.5) * h[i]);
      }

      local_sum += func(x);
    }
    return local_sum;
  }, std::plus<>());

  GetOutput() = sum * cell_v;

  return std::isfinite(GetOutput());
}

bool GalkinDMultidimIntegralsRectanglesTBB::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace galkin_d_multidim_integrals_rectangles
