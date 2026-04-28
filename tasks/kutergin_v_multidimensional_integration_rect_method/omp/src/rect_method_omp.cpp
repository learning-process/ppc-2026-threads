#include "../include/rect_method_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "../../common/include/common.hpp"
#include "util/include/util.hpp"

namespace kutergin_v_multidimensional_integration_rect_method {

RectMethodOMP::RectMethodOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool RectMethodOMP::ValidationImpl() {
  const auto &input = GetInput();
  if (input.limits.size() != input.n_steps.size() || input.limits.empty()) {
    return false;
  }
  return std::ranges::all_of(input.n_steps, [](int n) { return n > 0; });
}

bool RectMethodOMP::PreProcessingImpl() {
  local_input_ = GetInput();
  res_ = 0.0;
  return true;
}

bool RectMethodOMP::RunImpl() {
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

  double total_sum = 0.0;

  omp_set_num_threads(ppc::util::GetNumThreads());  // установка числа потоков

  // параллельный регион
#pragma omp parallel default(none) shared(total_iterations, h) \
    reduction(+ : total_sum)  // reduction(+:total_sum) автоматически соберет сумму со всех потоков
  {
    size_t tid = omp_get_thread_num();    // индекс потока
    int threads = omp_get_num_threads();  // число потоков

    size_t chunk = total_iterations / threads;
    size_t remainder = total_iterations % threads;

    // вычисление стартового индекса и количество работы для текущего потока
    size_t my_start = (tid * chunk) + std::min(tid, remainder);
    size_t my_count = chunk + (tid < remainder ? 1 : 0);

    total_sum += CalculateChunkSum(my_start, my_count, h);
  }

  res_ = total_sum * d_v;
  return true;
}

bool RectMethodOMP::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

double RectMethodOMP::CalculateChunkSum(size_t start_idx, size_t count, const std::vector<double> &h) {
  if (count == 0) {
    return 0.0;
  }

  size_t dims = local_input_.limits.size();  // число размерностей пространства
  std::vector<int> current_indices(dims, 0);
  std::vector<double> coords(dims);  // создание вектора координат размером dims

  size_t temp_idx = start_idx;  // временный линейный индекс
  for (int dim = static_cast<int>(dims) - 1; dim >= 0; --dim) {
    current_indices[dim] = static_cast<int>(temp_idx % local_input_.n_steps[dim]);  // текущий индекс в измерении d
    temp_idx /= local_input_.n_steps[dim];                                          // переход к следующему измерению
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
