#include "savva_d_monte_carlo/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "savva_d_monte_carlo/common/include/common.hpp"

namespace savva_d_monte_carlo {

SavvaDMonteCarloTBB::SavvaDMonteCarloTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SavvaDMonteCarloTBB::ValidationImpl() {
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

bool SavvaDMonteCarloTBB::PreProcessingImpl() {
  return true;
}

bool SavvaDMonteCarloTBB::RunImpl() {
  const auto &input = GetInput();
  auto &result = GetOutput();

  const size_t dim = input.Dimension();
  const double vol = input.Volume();
  const auto n = static_cast<int64_t>(input.count_points);
  const auto &func = input.f;

  const auto &lb = input.lower_bounds;
  const auto &ub = input.upper_bounds;

  std::vector<std::uniform_real_distribution<double>> dists;
  dists.reserve(dim);
  for (size_t i = 0; i < dim; ++i) {
    dists.emplace_back(lb[i], ub[i]);
  }

  double sum = tbb::parallel_reduce(tbb::blocked_range<int64_t>(0, n), 0.0,

                                    [&](const tbb::blocked_range<int64_t> &r, double local_sum) {
    std::minstd_rand gen(static_cast<uint32_t>(r.begin()) ^ 0x9e3779b9U);

    std::vector<double> point(dim);

    for (int64_t i = r.begin(); i < r.end(); ++i) {
      // генерация точки
      for (size_t j = 0; j < dim; ++j) {
        point[j] = dists[j](gen);
      }

      local_sum += func(point);
    }

    return local_sum;
  },

                                    [](double a, double b) { return a + b; });

  result = vol * sum / static_cast<double>(n);
  return true;
}

bool SavvaDMonteCarloTBB::PostProcessingImpl() {
  return true;
}

}  // namespace savva_d_monte_carlo
