#pragma once

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace nikitin_a_monte_carlo {

// Перечисление для выбора тестовой функции
enum class FunctionType : std::uint8_t {
  kConstant,
  kLinear,
  kProduct,
  kQuadratic,
  kExponential,
};

// Входные данные: (нижние границы, верхние границы, количество точек, тип функции)
using InType = std::tuple<std::vector<double>, std::vector<double>, int, FunctionType>;

// Выходные данные: приближенное значение интеграла
using OutType = double;

// Для тестов: (входные данные, ожидаемый результат в виде строки)
using TestType = std::tuple<int, std::string, InType>;

// Базовый класс задачи
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace nikitin_a_monte_carlo
