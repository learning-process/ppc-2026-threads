#include "nikitin_a_monte_carlo/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

#include "nikitin_a_monte_carlo/common/include/common.hpp"

namespace nikitin_a_monte_carlo {

namespace {
double EvaluateFunction(const std::vector<double> &point, FunctionType type) {
  if (point.empty()) {
    return 0.0;
  }
  switch (type) {
    case FunctionType::kConstant:
      return 1.0;
    case FunctionType::kLinear:
      return point[0];
    case FunctionType::kProduct:
      if (point.size() >= 2) {
        return point[0] * point[1];
      }
      return 0.0;
    case FunctionType::kQuadratic:
      if (point.size() >= 2) {
        return (point[0] * point[0]) + (point[1] * point[1]);
      }
      return 0.0;
    case FunctionType::kExponential:
      return std::exp(point[0]);
    default:
      return 0.0;
  }
}

double KroneckerSequence(int index, int dimension) {
  const std::array<double, 10> primes = {2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 19.0, 23.0, 29.0};
  const auto dim_idx = static_cast<std::size_t>(dimension % 10);
  // Используем at() вместо [] для проверки границ
  double alpha = std::sqrt(primes.at(dim_idx));
  alpha = alpha - std::floor(alpha);
  return std::fmod(static_cast<double>(index) * alpha, 1.0);
}
}  // namespace

NikitinAMonteCarloALL::NikitinAMonteCarloALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool NikitinAMonteCarloALL::ValidationImpl() {
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

bool NikitinAMonteCarloALL::PreProcessingImpl() {
  return true;
}

bool NikitinAMonteCarloALL::PostProcessingImpl() {
  return true;
}

bool NikitinAMonteCarloALL::RunImpl() {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const auto input = GetInput();
  const auto &lower_bounds = std::get<0>(input);
  const auto &upper_bounds = std::get<1>(input);
  const int num_points = std::get<2>(input);
  const FunctionType func_type = std::get<3>(input);
  const std::size_t dim = lower_bounds.size();

  double volume = 1.0;
  for (std::size_t i = 0; i < dim; ++i) {
    volume *= (upper_bounds[i] - lower_bounds[i]);
  }

  int points_per_rank = num_points / size;
  int remainder = num_points % size;

  int my_start = (rank * points_per_rank) + std::min(rank, remainder);
  int my_end = my_start + points_per_rank + (rank < remainder ? 1 : 0);

  double local_sum = 0.0;

#pragma omp parallel for reduction(+ : local_sum) schedule(static) default(none) \
    shared(my_start, my_end, dim, lower_bounds, upper_bounds, func_type)
  for (int i = my_start; i < my_end; ++i) {
    std::vector<double> point(dim);
    for (std::size_t j = 0; j < dim; ++j) {
      double u = KroneckerSequence(i, static_cast<int>(j));
      point[j] = lower_bounds[j] + (u * (upper_bounds[j] - lower_bounds[j]));
    }
    local_sum += EvaluateFunction(point, func_type);
  }

  double global_sum = 0.0;
  MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    double result = volume * global_sum / static_cast<double>(num_points);
    GetOutput() = result;
  }

  MPI_Bcast(&GetOutput(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  std::atomic<int> counter(0);
  std::thread t([&]() { counter++; });
  t.join();

  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

}  // namespace nikitin_a_monte_carlo
