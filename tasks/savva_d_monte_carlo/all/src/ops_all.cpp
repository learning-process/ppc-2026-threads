#include "savva_d_monte_carlo/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <cstddef>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

#include "savva_d_monte_carlo/common/include/common.hpp"

namespace savva_d_monte_carlo {

SavvaDMonteCarloALL::SavvaDMonteCarloALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SavvaDMonteCarloALL::ValidationImpl() {
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

bool SavvaDMonteCarloALL::PreProcessingImpl() {
  return true;
}

bool SavvaDMonteCarloALL::RunImpl() {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  const auto &input = GetInput();

  const size_t dim = input.Dimension();
  const uint64_t n_total = input.count_points;

  // 1. Распределение нагрузки
  uint64_t n_local = n_total / size;
  uint64_t remainder = n_total % size;
  if (std::cmp_less(rank, remainder)) {
    n_local++;
  }

  double local_sum = 0.0;
  const auto &func = input.f;

// 2. Параллельные вычисления внутри узла через OMP
#pragma omp parallel default(none) shared(input, func, dim, n_local, rank) reduction(+ : local_sum)
  {
    auto seed = static_cast<uint32_t>(rank * omp_get_num_threads() + omp_get_thread_num());
    std::minstd_rand gen(seed ^ 0x9e3779b9U);

    std::vector<std::uniform_real_distribution<double>> dists(dim);
    for (size_t i = 0; i < dim; ++i) {
      dists[i] = std::uniform_real_distribution<double>(input.lower_bounds[i], input.upper_bounds[i]);
    }

    std::vector<double> point(dim);

#pragma omp for schedule(static)
    for (int64_t i = 0; i < n_local; ++i) {
      for (size_t j = 0; j < dim; ++j) {
        point[j] = dists[j](gen);
      }
      local_sum += func(point);
    }
  }

  // 3. Глобальная редукция через MPI
  double total_sum = 0.0;
  MPI_Reduce(&local_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    double vol = input.Volume();
    GetOutput() = (total_sum / static_cast<double>(n_total)) * vol;
  }
  MPI_Bcast(&GetOutput(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  return true;
}

bool SavvaDMonteCarloALL::PostProcessingImpl() {
  return true;
}

}  // namespace savva_d_monte_carlo
