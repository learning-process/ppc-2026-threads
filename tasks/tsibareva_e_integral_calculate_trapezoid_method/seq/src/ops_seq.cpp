#include "tsibareva_e_integral_calculate_trapezoid_method/seq/include/ops_seq.hpp"
#include "tsibareva_e_integral_calculate_trapezoid_method/common/include/common.hpp"

#include <cmath>
#include <functional>
#include <vector>
#include <stdexcept>

namespace tsibareva_e_integral_calculate_trapezoid_method {

TsibarevaEIntegralCalculateTrapezoidMethodSEQ::TsibarevaEIntegralCalculateTrapezoidMethodSEQ(const IntegralInput& in)
    : ppc::task::Task<IntegralInput, double>() {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodSEQ::ValidationImpl() {
  return true;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodSEQ::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodSEQ::RunImpl() {
  const auto& input = GetInput();
  int dim = input.dimension;

    // Вычисляем шаги по каждому измерению
    std::vector<double> h(dim);
    for (int d = 0; d < dim; ++d) {
      h[d] = (input.upper_bounds[d] - input.lower_bounds[d]) / input.num_steps[d];
    }
    
    // Индексы текущей точки
    std::vector<int> indices(dim, 0);
    double sum = 0.0;
    
    // Перебор всех узлов сетки
    while (true) {
      // Формируем точку
      std::vector<double> point(dim);
      for (int d = 0; d < dim; ++d) {
        point[d] = input.lower_bounds[d] + indices[d] * h[d];
      }
      
      // Вычисляем вес для метода трапеций
      int boundary_count = 0;
      for (int d = 0; d < dim; ++d) {
        if (indices[d] == 0 || indices[d] == input.num_steps[d]) {
          ++boundary_count;
        }
      }
      
      // Вес = 0.5 для каждой граничной координаты
      double weight = (boundary_count == 0) ? 1.0 : std::pow(0.5, boundary_count);
      sum += weight * input.function(point);
      
      // Переход к следующей точке
      int d = dim - 1;
      while (d >= 0) {
        if (++indices[d] <= input.num_steps[d]) {
          break;
        }
        indices[d] = 0;
        --d;
      }
      if (d < 0) break;
    }
    
    // Умножаем на произведение шагов
    double product_of_steps = 1.0;
    for (int d = 0; d < dim; ++d) {
      product_of_steps *= h[d];
    }
    
    GetOutput() = sum * product_of_steps;
  
  return true;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace tsibareva_e_integral_calculate_trapezoid_method
