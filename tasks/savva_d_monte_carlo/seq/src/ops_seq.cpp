#include "savva_d_monte_carlo/seq/include/ops_seq.hpp"

#include <cstdint>
#include <random>
#include <vector>

#include "savva_d_monte_carlo/common/include/common.hpp"

namespace savva_d_monte_carlo {

SavvaDMonteCarloSEQ::SavvaDMonteCarloSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SavvaDMonteCarloSEQ::ValidationImpl() {
  const auto &input = GetInput();

  // Проверка количества точек
  if (input.count_points == 0) {
    return false;
  }

  // Проверка наличия функции
  if (!input.f) {
    return false;
  }

  // Проверка размерности
  if (input.Dimension() == 0) {
    return false;
  }

  // Проверка корректности границ
  for (size_t i = 0; i < input.Dimension(); ++i) {
    if (input.lower_bounds[i] > input.upper_bounds[i]) {
      return false;
    }
  }

  return true;
}

bool SavvaDMonteCarloSEQ::PreProcessingImpl() {
  return true;
}

bool SavvaDMonteCarloSEQ::RunImpl() {
  const auto &input = GetInput();
  auto &result = GetOutput();
  static thread_local std::minstd_rand generator(std::random_device{}());

  const size_t dim = input.Dimension();
  const double vol = input.Volume();
  const uint64_t n = input.count_points;
  const auto &func = input.f;

  std::vector<std::uniform_real_distribution<double>> distributions;
  distributions.resize(dim);

  for (size_t i = 0; i < dim; ++i) {
    distributions[i] = std::uniform_real_distribution<double>(input.lower_bounds[i], input.upper_bounds[i]);
  }

  double sum = 0.0;
  uint64_t i = 0;

  std::vector<double> p1(dim);
  std::vector<double> p2(dim);
  std::vector<double> p3(dim);
  std::vector<double> p4(dim);

  for (; i + 3 < n; i += 4) {
    for (size_t j = 0; j < dim; ++j) {
      p1[j] = distributions[j](generator);
      p2[j] = distributions[j](generator);
      p3[j] = distributions[j](generator);
      p4[j] = distributions[j](generator);
    }
    sum += func(p1) + func(p2) + func(p3) + func(p4);
  }

  // Обрабатываем оставшиеся точки
  for (; i < n; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      p1[j] = distributions[j](generator);
    }
    sum += func(p1);
  }

  // Вычисляем среднее и умножаем на объем (численно устойчивый способ)
  double mean = sum / static_cast<double>(n);
  result = mean * vol;

  return true;
}

bool SavvaDMonteCarloSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace savva_d_monte_carlo
