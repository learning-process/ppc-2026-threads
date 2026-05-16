#include "galkin_d_multidim_integrals_rectangles/stl/include/ops_stl.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <thread>
#include <utility>
#include <vector>

#include "galkin_d_multidim_integrals_rectangles/common/include/common.hpp"
#include "util/include/util.hpp"

namespace {

using Borders = std::vector<std::pair<double, double>>;

bool ComputeGridParams(const Borders &borders, int n, std::vector<double> &h, double &cell_v) {
  for (std::size_t i = 0; i < borders.size(); ++i) {
    h[i] = (borders[i].second - borders[i].first) / static_cast<double>(n);
    if (!(h[i] > 0.0) || !std::isfinite(h[i])) {
      return false;
    }
    cell_v *= h[i];
  }
  return true;
}

bool ComputeTotalCells(std::size_t dim, int n, std::int64_t &out) {
  std::size_t total = 1;
  for (std::size_t i = 0; i < dim; ++i) {
    if (total > (std::numeric_limits<std::size_t>::max() / static_cast<std::size_t>(n))) {
      return false;
    }
    total *= static_cast<std::size_t>(n);
  }
  if (total > static_cast<std::size_t>(LLONG_MAX)) {
    return false;
  }
  out = static_cast<std::int64_t>(total);
  return true;
}

}  // namespace

namespace galkin_d_multidim_integrals_rectangles {

GalkinDMultidimIntegralsRectanglesSTL::GalkinDMultidimIntegralsRectanglesSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool GalkinDMultidimIntegralsRectanglesSTL::ValidationImpl() {
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

bool GalkinDMultidimIntegralsRectanglesSTL::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool GalkinDMultidimIntegralsRectanglesSTL::RunImpl() {
  const auto &[func, borders, n] = GetInput();
  const std::size_t dim = borders.size();

  std::vector<double> h(dim);
  double cell_v = 1.0;
  if (!ComputeGridParams(borders, n, h, cell_v)) {
    return false;
  }

  std::int64_t total_i64 = 0;
  if (!ComputeTotalCells(dim, n, total_i64)) {
    return false;
  }
  const int requested = ppc::util::GetNumThreads();
  const int positive = (requested > 0) ? requested : 1;
  const int cap = (total_i64 > static_cast<std::int64_t>(std::numeric_limits<int>::max()))
                      ? std::numeric_limits<int>::max()
                      : static_cast<int>(total_i64);
  const int num_threads = std::min(positive, cap);

  std::vector<double> partial_sums(num_threads, 0.0);
  const std::int64_t chunk = total_i64 / num_threads;
  const std::int64_t remainder = total_i64 % num_threads;

  auto worker = [&](int tid) {
    const std::int64_t begin =
        (static_cast<std::int64_t>(tid) * chunk) + std::min(static_cast<std::int64_t>(tid), remainder);
    const std::int64_t end = begin + chunk + (static_cast<std::int64_t>(tid) < remainder ? 1 : 0);
    std::vector<double> x(dim);
    double local_sum = 0.0;
    for (std::int64_t linear_idx = begin; linear_idx < end; ++linear_idx) {
      auto tmp = static_cast<std::size_t>(linear_idx);
      for (std::size_t i = 0; i < dim; ++i) {
        const std::size_t idx_i = tmp % static_cast<std::size_t>(n);
        tmp /= static_cast<std::size_t>(n);
        x[i] = borders[i].first + ((static_cast<double>(idx_i) + 0.5) * h[i]);
      }
      local_sum += func(x);
    }
    partial_sums[tid] = local_sum;
  };

  std::vector<std::thread> threads(num_threads - 1);
  for (int ti = 1; ti < num_threads; ++ti) {
    threads[ti - 1] = std::thread(worker, ti);
  }
  worker(0);

  double sum = partial_sums[0];
  for (int ti = 1; ti < num_threads; ++ti) {
    threads[ti - 1].join();
    sum += partial_sums[ti];
  }

  GetOutput() = sum * cell_v;
  return std::isfinite(GetOutput());
}

bool GalkinDMultidimIntegralsRectanglesSTL::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace galkin_d_multidim_integrals_rectangles
