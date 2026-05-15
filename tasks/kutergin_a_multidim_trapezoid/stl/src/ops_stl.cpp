#include "kutergin_a_multidim_trapezoid/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <thread>
#include <vector>

#include "kutergin_a_multidim_trapezoid/common/include/common.hpp"

namespace kutergin_a_multidim_trapezoid {

KuterginAMultidimTrapezoidSTL::KuterginAMultidimTrapezoidSTL(const InType &in) {
  this->SetTypeOfTask(GetStaticTypeOfTask());
  this->GetInput() = in;
  this->GetOutput() = 0.0;
}

bool KuterginAMultidimTrapezoidSTL::ValidationImpl() {
  const auto &input = GetInput();
  return std::get<0>(input) != nullptr && std::get<2>(input) > 0 && !std::get<1>(input).empty();
}

bool KuterginAMultidimTrapezoidSTL::PreProcessingImpl() {
  res_ = 0.0;
  return true;
}

bool KuterginAMultidimTrapezoidSTL::RunImpl() {
  const auto &input = GetInput();
  const auto &borders = std::get<1>(input);
  const int n = std::get<2>(input);
  const size_t dims = borders.size();

  size_t total_nodes = 1;
  std::vector<double> h(dims);
  double cell_volume = 1.0;

  for (size_t i = 0; i < dims; ++i) {
    total_nodes *= (static_cast<size_t>(n) + 1);
    h[i] = (borders[i].second - borders[i].first) / n;
    cell_volume *= h[i];
  }

  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  if (total_nodes < static_cast<size_t>(num_threads)) {
    num_threads = static_cast<unsigned int>(total_nodes);
  }

  std::vector<std::thread> threads;
  std::vector<double> partial_sums(num_threads, 0.0);

  const size_t chunk = total_nodes / num_threads;
  const size_t remainder = total_nodes % num_threads;

  for (unsigned int i = 0; i < num_threads; ++i) {
    const size_t start = (i * chunk) + std::min(static_cast<size_t>(i), remainder);
    const size_t end = start + chunk + (static_cast<size_t>(i) < remainder ? 1 : 0);

    if (start < end) {
      threads.emplace_back([this, start, end, &h, &partial_sums, i]() noexcept {
        try {
          partial_sums[i] = CalculateChunkSum(start, end, h);
        } catch (...) {
          (void)0;
        }
      });
    }
  }

  for (auto &t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }

  const double total_weight_sum = std::accumulate(partial_sums.begin(), partial_sums.end(), 0.0);
  res_ = total_weight_sum * cell_volume;

  return std::isfinite(res_);
}

bool KuterginAMultidimTrapezoidSTL::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

double KuterginAMultidimTrapezoidSTL::CalculateChunkSum(size_t start_idx, size_t end_idx,
                                                        const std::vector<double> &h) {
  const auto &input = GetInput();
  const auto &func = std::get<0>(input);
  const auto &borders = std::get<1>(input);
  const int n = std::get<2>(input);
  const size_t dims = borders.size();

  std::vector<double> coords(dims);
  double chunk_sum = 0.0;

  for (size_t i = start_idx; i < end_idx; ++i) {
    double weight = 1.0;
    size_t temp_idx = i;

    for (size_t dim_idx = 0; dim_idx < dims; ++dim_idx) {
      const int coord_idx = static_cast<int>(temp_idx % (static_cast<size_t>(n) + 1));
      temp_idx /= (static_cast<size_t>(n) + 1);

      coords[dim_idx] = borders[dim_idx].first + (static_cast<double>(coord_idx) * h[dim_idx]);

      if (coord_idx == 0 || coord_idx == n) {
        weight *= 0.5;
      }
    }
    chunk_sum += weight * func(coords);
  }
  return chunk_sum;
}

}  // namespace kutergin_a_multidim_trapezoid
