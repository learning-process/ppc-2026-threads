#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace savva_d_monte_carlo {

struct InputData {
  std::vector<double> lower_bounds;
  std::vector<double> upper_bounds;
  uint64_t count_points = 0;                             // Количество точек для выборки
  std::function<double(const std::vector<double> &)> f;  // Функция от вектора

  InputData() = default;

  // Конструктор для проверки совпадения размерностей
  InputData(std::vector<double> lowers, std::vector<double> uppers, uint64_t n,
            std::function<double(const std::vector<double> &)> func)
      : lower_bounds(std::move(lowers)), upper_bounds(std::move(uppers)), count_points(n), f(std::move(func)) {
    // Проверяем, что размерности совпадают
    if (lower_bounds.size() != upper_bounds.size()) {
      throw std::invalid_argument("Lower and upper bounds must have the same dimension");
    }
    if (n == 0) {
      throw std::invalid_argument("Number of points must be positive");
    }
  }

  // Возвращает размерность задачи
  size_t Dimension() const {
    return lower_bounds.size();
  }

  // Вычисляет объем области интегрирования
  double Volume() const {
    if (lower_bounds.empty()) {
      return 0.0;
    }

    double vol = 1.0;
    for (size_t i = 0; i < lower_bounds.size(); ++i) {
      vol *= (upper_bounds[i] - lower_bounds[i]);
    }
    return vol;
  }
};

using InType = InputData;
using OutType = double;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace savva_d_monte_carlo
