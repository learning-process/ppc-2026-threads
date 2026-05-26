#include "nazyrov_a_multidim_integral_rectangle/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

#include "nazyrov_a_multidim_integral_rectangle/common/include/common.hpp"
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_reduce.h"

namespace nazyrov_a_multidim_integral_rectangle {

NazyrovAMultidimIntegralRectangleTbb::NazyrovAMultidimIntegralRectangleTbb(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool NazyrovAMultidimIntegralRectangleTbb::ValidationImpl() {
  const auto &func = std::get<0>(GetInput());
  const auto &bounds = std::get<1>(GetInput());
  const int n = std::get<2>(GetInput());
  return func && n > 0 && !bounds.empty() && std::ranges::all_of(bounds, [](const auto &bd) {
    return std::isfinite(bd.first) && std::isfinite(bd.second) && bd.first < bd.second;
  });
}

bool NazyrovAMultidimIntegralRectangleTbb::PreProcessingImpl() {
  return true;
}

bool NazyrovAMultidimIntegralRectangleTbb::RunImpl() {
  const auto &func = std::get<0>(GetInput());
  const auto &bounds = std::get<1>(GetInput());
  const int n = std::get<2>(GetInput());

  const int dim = static_cast<int>(bounds.size());

  std::vector<double> h(static_cast<std::size_t>(dim));
  double cell_vol = 1.0;
  for (int i = 0; i < dim; ++i) {
    h[i] = (bounds[i].second - bounds[i].first) / n;
    cell_vol *= h[i];
  }

  std::int64_t total = 1;
  for (int i = 0; i < dim; ++i) {
    total *= n;
  }

  const double sum =
      oneapi::tbb::parallel_reduce(oneapi::tbb::blocked_range<std::int64_t>(0, total), 0.0,
                                   [&](const oneapi::tbb::blocked_range<std::int64_t> &range, double partial) {
    std::vector<double> point(static_cast<std::size_t>(dim));
    for (std::int64_t cell = range.begin(); cell < range.end(); ++cell) {
      std::int64_t tmp = cell;
      for (int i = dim - 1; i >= 0; --i) {
        const int ki = static_cast<int>(tmp % n);
        tmp /= n;
        const double coordinate = bounds[i].first + ((ki + 0.5) * h[i]);
        point[i] = coordinate;
      }
      partial += func(point);
    }
    return partial;
  }, std::plus<double>{});

  GetOutput() = sum * cell_vol;
  return true;
}

bool NazyrovAMultidimIntegralRectangleTbb::PostProcessingImpl() {
  return true;
}

}  // namespace nazyrov_a_multidim_integral_rectangle
