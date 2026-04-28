#include "nikitin_a_monte_carlo/omp/include/ops_omp.hpp"

#include <omp.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "nikitin_a_monte_carlo/common/include/common.hpp"

namespace nikitin_a_monte_carlo {

namespace {
// Вспомогательная функция для вычисления значения тестовой функции
double EvaluateFunction(const std::vector<double> &point, FunctionType type) {
  if (point.empty()) {
    return 0.0;
  }

  switch (type) {
    case FunctionType::kConstant:
      return 1.0;
    case FunctionType::kLinear:
      return point.at(0);
    case FunctionType::kProduct:
      if (point.size() < 2) {
        return 0.0;
      }
      return point.at(0) * point.at(1);
    case FunctionType::kQuadratic:
      if (point.size() < 2) {
        return 0.0;
      }
      return (point.at(0) * point.at(0)) + (point.at(1) * point.at(1));
    case FunctionType::kExponential:
      return std::exp(point.at(0));
    default:
      return 0.0;
  }
}

// Генерация квазислучайной последовательности Кронекера
double KroneckerSequence(int index, int dimension) {
  const std::array<double, 10> primes = {2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 19.0, 23.0, 29.0};
  double alpha = std::sqrt(primes.at(static_cast<std::size_t>(dimension % 10)));
  alpha = alpha - std::floor(alpha);
  return std::fmod(index * alpha, 1.0);
}
}  // namespace

NikitinAMonteCarloOMP::NikitinAMonteCarloOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool NikitinAMonteCarloOMP::ValidationImpl() {
  const auto &[lower_bounds, upper_bounds, num_points, func_type] = GetInput();

  if (lower_bounds.empty() || upper_bounds.empty()) {
    return false;
  }

  if (lower_bounds.size() != upper_bounds.size()) {
    return false;
  }

  for (std::size_t i = 0; i < lower_bounds.size(); ++i) {
    if (lower_bounds[i] >= upper_bounds[i]) {
      return false;
    }
  }

  return num_points > 0;
}

bool NikitinAMonteCarloOMP::PreProcessingImpl() {
  return true;
}

bool NikitinAMonteCarloOMP::RunImpl() {
  // Получаем входные данные и распаковываем их в обычные переменные
  const auto input = GetInput();
  const auto &lower_bounds = std::get<0>(input);
  const auto &upper_bounds = std::get<1>(input);
  const int num_points = std::get<2>(input);
  const FunctionType func_type = std::get<3>(input);

  std::size_t dim = lower_bounds.size();

  // Вычисление объема области интегрирования
  double volume = 1.0;
  for (std::size_t i = 0; i < dim; ++i) {
    volume *= (upper_bounds[i] - lower_bounds[i]);
  }

  double sum = 0.0;

// Параллельный цикл с явным указанием всех используемых переменных
#pragma omp parallel for schedule(static) reduction(+ : sum) default(none) \
    shared(lower_bounds, upper_bounds, dim, num_points, func_type)
  for (int i = 0; i < num_points; ++i) {
    std::vector<double> point(dim);
    for (std::size_t j = 0; j < dim; ++j) {
      double u = KroneckerSequence(i, static_cast<int>(j));
      point[j] = lower_bounds[j] + (u * (upper_bounds[j] - lower_bounds[j]));
    }
    sum += EvaluateFunction(point, func_type);
  }

  double result = volume * sum / static_cast<double>(num_points);
  GetOutput() = result;

  return true;
}

bool NikitinAMonteCarloOMP::PostProcessingImpl() {
  return true;
}

}  // namespace nikitin_a_monte_carlo
