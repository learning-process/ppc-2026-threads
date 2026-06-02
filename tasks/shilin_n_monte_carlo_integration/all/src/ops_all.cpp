#include "shilin_n_monte_carlo_integration/all/include/ops_all.hpp"

#include <mpi.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "shilin_n_monte_carlo_integration/common/include/common.hpp"

namespace shilin_n_monte_carlo_integration {

ShilinNMonteCarloIntegrationALL::ShilinNMonteCarloIntegrationALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ShilinNMonteCarloIntegrationALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) {
    return true;
  }
  const auto &[lower, upper, n, func_type] = GetInput();
  if (lower.size() != upper.size() || lower.empty()) {
    return false;
  }
  if (n <= 0) {
    return false;
  }
  for (size_t i = 0; i < lower.size(); ++i) {
    if (lower[i] >= upper[i]) {
      return false;
    }
  }
  if (func_type < FuncType::kConstant || func_type > FuncType::kSinProduct) {
    return false;
  }
  constexpr size_t kMaxDimensions = 10;
  return lower.size() <= kMaxDimensions;
}

bool ShilinNMonteCarloIntegrationALL::PreProcessingImpl() {
  const auto &[lower, upper, n, func_type] = GetInput();
  lower_bounds_ = lower;
  upper_bounds_ = upper;
  num_points_ = n;
  func_type_ = func_type;
  return true;
}

bool ShilinNMonteCarloIntegrationALL::RunImpl() {
  int rank = 0;
  int num_ranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  auto dimensions = static_cast<int>(lower_bounds_.size());

  const std::vector<double> alpha = {
      0.41421356237309504,  // frac(sqrt(2))
      0.73205080756887729,  // frac(sqrt(3))
      0.23606797749978969,  // frac(sqrt(5))
      0.64575131106459059,  // frac(sqrt(7))
      0.31662479035539984,  // frac(sqrt(11))
      0.60555127546398929,  // frac(sqrt(13))
      0.12310562561766059,  // frac(sqrt(17))
      0.35889894354067355,  // frac(sqrt(19))
      0.79583152331271838,  // frac(sqrt(23))
      0.38516480713450403   // frac(sqrt(29))
  };

  int local_count = 0;
  if (rank < num_points_) {
    local_count = (((num_points_ - rank - 1) / num_ranks) + 1);
  }

  double local_sum = 0.0;

  // MSVC OpenMP does not allow non-static data members in data-sharing clauses; use `this`.
  auto *self = this;
#pragma omp parallel default(none) shared(dimensions, alpha, self, rank, num_ranks, local_count) \
    reduction(+ : local_sum)
  {
    std::vector<double> point(dimensions);
#pragma omp for schedule(static)
    for (int k = 0; k < local_count; ++k) {
      int i = rank + (k * num_ranks);
      for (int di = 0; di < dimensions; ++di) {
        double val = 0.5 + (static_cast<double>(i + 1) * alpha[di]);
        double current = val - std::floor(val);
        point[di] = self->lower_bounds_[di] + ((self->upper_bounds_[di] - self->lower_bounds_[di]) * current);
      }
      local_sum += IntegrandFunction::Evaluate(self->func_type_, point);
    }
  }

  double global_sum = 0.0;
  MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double volume = 1.0;
  for (int di = 0; di < dimensions; ++di) {
    volume *= (upper_bounds_[di] - lower_bounds_[di]);
  }

  GetOutput() = volume * global_sum / static_cast<double>(num_points_);
  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

bool ShilinNMonteCarloIntegrationALL::PostProcessingImpl() {
  return true;
}

}  // namespace shilin_n_monte_carlo_integration
