#include "nikitin_a_monte_carlo/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

#include "nikitin_a_monte_carlo/common/include/common.hpp"

namespace nikitin_a_monte_carlo {

namespace {
// Вспомогательная функция для вычисления значения тестовой функции
double EvaluateFunction(const std::vector<double> &point, FunctionType type) {
  switch (type) {
    case FunctionType::kConstant:
      return 1.0;
    case FunctionType::kLinear:
      return point[0];
    case FunctionType::kProduct:
      return point[0] * point[1];
    case FunctionType::kQuadratic:
      return point[0] * point[0] + point[1] * point[1];
    case FunctionType::kExponential:
      return std::exp(point[0]);
    default:
      return 0.0;
  }
}

// Генерация квазислучайной последовательности Кронекера
double KroneckerSequence(int index, int dimension) {
  // Используем простые числа для разных измерений
  const double primes[] = {2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 19.0, 23.0, 29.0};
  double alpha = std::sqrt(primes[dimension % 10]);
  // Берем дробную часть
  alpha = alpha - std::floor(alpha);
  return std::fmod(index * alpha, 1.0);
}
}  // namespace

NikitinAMonteCarloSEQ::NikitinAMonteCarloSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool NikitinAMonteCarloSEQ::ValidationImpl() {
  const auto &[lower_bounds, upper_bounds, num_points, func_type] = GetInput();

  // Проверка на пустые векторы
  if (lower_bounds.empty() || upper_bounds.empty()) {
    return false;
  }

  // Размерности нижних и верхних границ должны совпадать
  if (lower_bounds.size() != upper_bounds.size()) {
    return false;
  }

  // Проверка, что нижняя граница меньше верхней для каждого измерения
  for (std::size_t i = 0; i < lower_bounds.size(); ++i) {
    if (lower_bounds[i] >= upper_bounds[i]) {
      return false;
    }
  }

  // Количество точек должно быть положительным
  if (num_points <= 0) {
    return false;
  }

  return true;
}

bool NikitinAMonteCarloSEQ::PreProcessingImpl() {
  // Предобработка не требуется, так как все данные уже в GetInput()
  return true;
}

bool NikitinAMonteCarloSEQ::RunImpl() {
  const auto &[lower_bounds, upper_bounds, num_points, func_type] = GetInput();

  std::size_t dim = lower_bounds.size();

  // Вычисление объема области интегрирования
  double volume = 1.0;
  for (std::size_t i = 0; i < dim; ++i) {
    volume *= (upper_bounds[i] - lower_bounds[i]);
  }

  // Сумма значений функции в точках
  double sum = 0.0;

  // Генерация точек и вычисление функции
  for (int i = 0; i < num_points; ++i) {
    std::vector<double> point(dim);

    // Генерация точки с помощью последовательности Кронекера
    for (std::size_t j = 0; j < dim; ++j) {
      // Получаем значение в единичном гиперкубе [0,1)
      double u = KroneckerSequence(i, static_cast<int>(j));
      // Масштабируем в реальную область интегрирования
      point[j] = lower_bounds[j] + u * (upper_bounds[j] - lower_bounds[j]);
    }

    // Вычисляем значение функции в точке
    sum += EvaluateFunction(point, func_type);
  }

  // Вычисляем приближенное значение интеграла
  double result = volume * sum / static_cast<double>(num_points);

  // Сохраняем результат
  GetOutput() = result;

  return true;
}

bool NikitinAMonteCarloSEQ::PostProcessingImpl() {
  // Постобработка не требуется, результат уже сохранен в GetOutput()
  return true;
}

}  // namespace nikitin_a_monte_carlo
