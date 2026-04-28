#include "../include/rect_method_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <thread>
#include <utility>
#include <vector>

#include "../../common/include/common.hpp"
#include "util/include/util.hpp"

namespace kutergin_v_multidimensional_integration_rect_method {

RectMethodSTL::RectMethodSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool RectMethodSTL::ValidationImpl() {
  const auto &input = GetInput();
  if (input.limits.size() != input.n_steps.size() || input.limits.empty()) {
    return false;
  }
  return std::ranges::all_of(input.n_steps, [](int n) { return n > 0; });
}

bool RectMethodSTL::PreProcessingImpl() {
  local_input_ = GetInput();
  res_ = 0.0;
  return true;
}

bool RectMethodSTL::RunImpl() {
  size_t dims = local_input_.limits.size();  // число размерностей пространства

  // вычисление числа шагов и общего числа итераций
  size_t total_iterations = 1;
  std::vector<double> h(dims);
  double d_v = 1.0;

  for (size_t i = 0; i < dims; i++) {
    total_iterations *= local_input_.n_steps[i];
    h[i] = (local_input_.limits[i].second - local_input_.limits[i].first) / local_input_.n_steps[i];
    d_v *= h[i];
  }

  int num_threads = ppc::util::GetNumThreads();

  std::vector<std::thread> threads;
  std::vector<double> partial_sums(num_threads, 0.0);

  size_t chunk = total_iterations / num_threads;
  size_t remainder = total_iterations % num_threads;

  for (int i = 0; i < num_threads; ++i) {
    size_t start = (i * chunk) + std::min(static_cast<size_t>(i), remainder);
    size_t end = start + chunk + (std::cmp_less(static_cast<size_t>(i), remainder) ? 1 : 0);

    if (start < end) {
      threads.emplace_back(
          [this, start, end, &h, &partial_sums, i]() { partial_sums[i] = CalculateChunkSum(start, end, h); });
    }
  }

  for (auto &t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }

  double total_sum = std::accumulate(partial_sums.begin(), partial_sums.end(), 0.0);
  res_ = total_sum * d_v;
  return true;
}

bool RectMethodSTL::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

double RectMethodSTL::CalculateChunkSum(size_t start_idx, size_t end_idx, const std::vector<double> &h) {
  if (start_idx >= end_idx) {
    return 0.0;
  }

  size_t count = end_idx - start_idx;
  size_t dims = local_input_.limits.size();  // число размерностей пространства
  std::vector<int> current_indices(dims, 0);
  std::vector<double> coords(dims);  // создание вектора координат размером dims

  size_t temp_idx = start_idx;
  for (int dim = static_cast<int>(dims) - 1; dim >= 0; --dim) {
    current_indices[dim] = static_cast<int>(temp_idx % local_input_.n_steps[dim]);
    temp_idx /= local_input_.n_steps[dim];
  }

  double chunk_sum = 0.0;

  for (size_t i = 0; i < count; ++i) {
    for (size_t dm = 0; dm < dims; ++dm) {
      coords[dm] = local_input_.limits[dm].first + ((current_indices[dm] + 0.5) * h[dm]);  // реальная координата
    }

    chunk_sum += local_input_.func(coords);  // вычисление функции в точке

    for (int dm = static_cast<int>(dims) - 1; dm >= 0; --dm) {
      current_indices[dm]++;
      if (current_indices[dm] < local_input_.n_steps[dm]) {
        break;
      }
      current_indices[dm] = 0;
    }
  }

  return chunk_sum;
}

}  // namespace kutergin_v_multidimensional_integration_rect_method
