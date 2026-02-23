#include "../include/rect_method_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "../../common/include/common.hpp"

namespace kutergin_v_multidimensional_integration_rect_method {

RectMethodSequential::RectMethodSequential(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool RectMethodSequential::ValidationImpl() {
  const auto &input = GetInput();
  if (input.limits.size() != input.n_steps.size() || input.limits.empty()) {
    return false;
  }
  return std::ranges::all_of(input.n_steps, [](int n) { return n > 0; });
}

bool RectMethodSequential::PreProcessingImpl() {
  local_input_ = GetInput();
  res_ = 0.0;
  return true;
}

bool RectMethodSequential::RunImpl() {
  size_t dims = local_input_.limits.size();  // число размерностей пространства
  std::vector<double> coords(dims);          // создание вектора координат размером dims

  size_t total_iterations = 1;
  std::vector<double> h(dims);

  double d_v = 1.0;
  for (size_t i = 0; i < dims; i++) {
    total_iterations *= local_input_.n_steps[i];
    h[i] = (local_input_.limits[i].second - local_input_.limits[i].first) / local_input_.n_steps[i];
    d_v *= h[i];
  }

  double total_sum = 0.0;
  std::vector<int> current_indices(dims, 0);

  for (size_t i = 0; i < total_iterations; ++i) {
    for (size_t dm = 0; dm < dims; ++dm) {
      coords[dm] = local_input_.limits[dm].first + ((current_indices[dm] + 0.5) * h[dm]);
    }
    total_sum += local_input_.func(coords);
    for (int dm = static_cast<int>(dims) - 1; dm >= 0; --dm) {
      current_indices[dm]++;
      if (current_indices[dm] < local_input_.n_steps[dm]) {
        break;
      }
      current_indices[dm] = 0;
    }
  }

  res_ = total_sum * d_v;
  return true;
}

bool RectMethodSequential::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

}  // namespace kutergin_v_multidimensional_integration_rect_method
