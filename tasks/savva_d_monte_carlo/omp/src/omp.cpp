#include <omp.h>

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "savva_d_monte_carlo/common/include/common.hpp"
#include "savva_d_monte_carlo/omp/include/ops_omp.hpp"

namespace savva_d_monte_carlo {

SavvaDMonteCarloOMP::SavvaDMonteCarloOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SavvaDMonteCarloOMP::ValidationImpl() {
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

bool SavvaDMonteCarloOMP::PreProcessingImpl() {
  return true;
}

bool SavvaDMonteCarloOMP::RunImpl() {
  const auto &input = GetInput();
  auto &result = GetOutput();

  const size_t dim = input.Dimension();
  const double vol = input.Volume();
  const auto n = static_cast<int64_t>(input.count_points);
  const auto &func = input.f;

  double sum = 0.0;

#pragma omp parallel default(none) shared(input, func, dim, n) reduction(+ : sum)
  {
    std::minstd_rand gen(1337 + omp_get_thread_num());

    std::vector<std::uniform_real_distribution<double>> dists(dim);
    for (size_t i = 0; i < dim; ++i) {
      dists[i] = std::uniform_real_distribution<double>(input.lower_bounds[i], input.upper_bounds[i]);
    }

    std::vector<double> point(dim);

#pragma omp for schedule(static)
    for (int64_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < dim; ++j) {
        point[j] = dists[j](gen);
      }
      sum += func(point);
    }
  }

  result = vol * sum / static_cast<double>(n);
  return true;
}

bool SavvaDMonteCarloOMP::PostProcessingImpl() {
  return true;
}

}  // namespace savva_d_monte_carlo
