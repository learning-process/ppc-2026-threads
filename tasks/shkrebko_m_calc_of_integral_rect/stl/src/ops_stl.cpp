#include "shkrebko_m_calc_of_integral_rect/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include "shkrebko_m_calc_of_integral_rect/common/include/common.hpp"
#include "util/include/util.hpp"

namespace shkrebko_m_calc_of_integral_rect {

ShkrebkoMCalcOfIntegralRectSTL::ShkrebkoMCalcOfIntegralRectSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ShkrebkoMCalcOfIntegralRectSTL::ValidationImpl() {
  const auto &input = GetInput();

  if (!input.func) {
    return false;
  }
  if (input.limits.size() != input.n_steps.size() || input.limits.empty()) {
    return false;
  }
  if (!std::ranges::all_of(input.n_steps, [](int n) { return n > 0; })) {
    return false;
  }
  if (!std::ranges::all_of(input.limits, [](const auto &lim) { return lim.first < lim.second; })) {
    return false;
  }
  return true;
}

bool ShkrebkoMCalcOfIntegralRectSTL::PreProcessingImpl() {
  local_input_ = GetInput();
  res_ = 0.0;
  return true;
}

namespace {

struct ChunkData {
  std::size_t start;
  std::size_t end;
  const std::vector<double> *h;
  const InType *input;
  double *output;
};

void ComputeChunkSum(ChunkData data) {
  const std::size_t dim = data.input->limits.size();
  const auto &limits = data.input->limits;
  const auto &n_steps = data.input->n_steps;
  const auto &func = data.input->func;
  const std::vector<double> &h = *data.h;

  thread_local std::vector<double> point;
  if (point.size() != dim) {
    point.resize(dim);
  }

  double local_sum = 0.0;
  for (std::size_t idx = data.start; idx < data.end; ++idx) {
    std::size_t tmp = idx;
    for (int i = static_cast<int>(dim) - 1; i >= 0; --i) {
      std::size_t coord_index = tmp % static_cast<std::size_t>(n_steps[i]);
      tmp /= static_cast<std::size_t>(n_steps[i]);
      point[i] = limits[i].first + ((static_cast<double>(coord_index) + 0.5) * h[i]);
    }
    local_sum += func(point);
  }
  *data.output = local_sum;
}

}  // namespace

bool ShkrebkoMCalcOfIntegralRectSTL::RunImpl() {
  const std::size_t dim = local_input_.limits.size();
  const auto &limits = local_input_.limits;
  const auto &n_steps = local_input_.n_steps;

  std::vector<double> h(dim);
  double cell_volume = 1.0;
  std::size_t total_points = 1;
  for (std::size_t i = 0; i < dim; ++i) {
    double left = limits[i].first;
    double right = limits[i].second;
    int steps = n_steps[i];
    h[i] = (right - left) / static_cast<double>(steps);
    cell_volume *= h[i];
    total_points *= static_cast<std::size_t>(steps);
  }

  int thread_count = ppc::util::GetNumThreads();
  thread_count = std::max(1, std::min(thread_count, static_cast<int>(total_points)));

  std::size_t chunk = total_points / thread_count;
  std::size_t remainder = total_points % thread_count;
  std::vector<std::thread> threads;
  std::vector<double> partial_sums(thread_count, 0.0);

  std::size_t start = 0;
  for (int i = 0; i < thread_count; ++i) {
    std::size_t extra = std::cmp_less(static_cast<std::size_t>(i), remainder) ? 1 : 0;
    std::size_t end = start + chunk + extra;
    if (start >= end) {
      partial_sums[i] = 0.0;
      continue;
    }
    ChunkData data = {.start = start, .end = end, .h = &h, .input = &local_input_, .output = &partial_sums[i]};
    threads.emplace_back(ComputeChunkSum, data);
    start = end;
  }

  for (auto &t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }

  double total_sum = 0.0;
  for (double s : partial_sums) {
    total_sum += s;
  }

  res_ = total_sum * cell_volume;
  return true;
}

bool ShkrebkoMCalcOfIntegralRectSTL::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

}  // namespace shkrebko_m_calc_of_integral_rect
