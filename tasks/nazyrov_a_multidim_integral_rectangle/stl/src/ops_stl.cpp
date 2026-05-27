#include "nazyrov_a_multidim_integral_rectangle/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <vector>

#include "nazyrov_a_multidim_integral_rectangle/common/include/common.hpp"
#include "util/include/util.hpp"

namespace nazyrov_a_multidim_integral_rectangle {

NazyrovAMultidimIntegralRectangleStl::NazyrovAMultidimIntegralRectangleStl(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool NazyrovAMultidimIntegralRectangleStl::ValidationImpl() {
  const auto &func = std::get<0>(GetInput());
  const auto &bounds = std::get<1>(GetInput());
  const int n = std::get<2>(GetInput());
  return func && n > 0 && !bounds.empty() && std::ranges::all_of(bounds, [](const auto &bd) {
    return std::isfinite(bd.first) && std::isfinite(bd.second) && bd.first < bd.second;
  });
}

bool NazyrovAMultidimIntegralRectangleStl::PreProcessingImpl() {
  return true;
}

bool NazyrovAMultidimIntegralRectangleStl::RunImpl() {
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

  const int num_threads = std::max(1, ppc::util::GetNumThreads());
  const int actual_threads = std::max(1, std::min(num_threads, static_cast<int>(total)));

  std::vector<double> partial_sums(static_cast<std::size_t>(actual_threads), 0.0);

  auto worker = [&](int tid, std::int64_t begin, std::int64_t end) {
    std::vector<double> point(static_cast<std::size_t>(dim));
    double local_sum = 0.0;
    for (std::int64_t cell = begin; cell < end; ++cell) {
      std::int64_t tmp = cell;
      for (int i = dim - 1; i >= 0; --i) {
        const int ki = static_cast<int>(tmp % n);
        tmp /= n;
        const double coordinate = bounds[i].first + ((ki + 0.5) * h[i]);
        point[i] = coordinate;
      }
      local_sum += func(point);
    }
    partial_sums[tid] = local_sum;
  };

  const std::int64_t chunk = total / actual_threads;
  const std::int64_t leftover = total % actual_threads;

  std::vector<std::thread> threads;
  threads.reserve(static_cast<std::size_t>(actual_threads));

  std::int64_t cursor = 0;
  for (int thread_idx = 0; thread_idx < actual_threads; ++thread_idx) {
    const std::int64_t begin = cursor;
    const std::int64_t end = begin + chunk + (static_cast<std::int64_t>(thread_idx) < leftover ? 1LL : 0LL);
    cursor = end;
    threads.emplace_back(worker, thread_idx, begin, end);
  }
  for (auto &th : threads) {
    th.join();
  }

  double sum = 0.0;
  for (double ps : partial_sums) {
    sum += ps;
  }

  GetOutput() = sum * cell_vol;
  return true;
}

bool NazyrovAMultidimIntegralRectangleStl::PostProcessingImpl() {
  return true;
}

}  // namespace nazyrov_a_multidim_integral_rectangle
