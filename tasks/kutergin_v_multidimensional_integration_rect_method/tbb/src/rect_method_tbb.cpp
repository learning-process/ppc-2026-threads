#include "../include/rect_method_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "../../common/include/common.hpp"
#include "util/include/util.hpp"

namespace kutergin_v_multidimensional_integration_rect_method {

RectMethodTBB::RectMethodTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool RectMethodTBB::ValidationImpl() {
  const auto &input = GetInput();
  if (input.limits.size() != input.n_steps.size() || input.limits.empty()) {
    return false;
  }
  return std::ranges::all_of(input.n_steps, [](int n) { return n > 0; });
}

bool RectMethodTBB::PreProcessingImpl() {
  local_input_ = GetInput();
  res_ = 0.0;
  return true;
}

bool RectMethodTBB::RunImpl() {
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
  tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, num_threads);

  double total_sum = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, total_iterations), 0.0,
                                          [&](const tbb::blocked_range<size_t> &r, double sum) {
    return sum + CalculateChunkSum(r.begin(), r.end(), h);
  }, std::plus<>());

  res_ = total_sum * d_v;
  return true;
}

bool RectMethodTBB::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

double RectMethodTBB::CalculateChunkSum(size_t start_idx, size_t end_idx, const std::vector<double> &h) {
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
