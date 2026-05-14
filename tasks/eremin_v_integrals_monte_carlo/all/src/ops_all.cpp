#include "eremin_v_integrals_monte_carlo/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

#include "eremin_v_integrals_monte_carlo/common/include/common.hpp"

namespace eremin_v_integrals_monte_carlo {

EreminVIntegralsMonteCarloALL::EreminVIntegralsMonteCarloALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool EreminVIntegralsMonteCarloALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int valid = 1;
  if (rank == 0) {
    const auto &input = GetInput();
    const bool bad = (input.samples <= 0) || input.bounds.empty() || (input.func == nullptr) ||
                     !std::ranges::all_of(input.bounds, [](const auto &p) {
      const auto &[a, b] = p;
      return (a < b) && (std::abs(a) <= 1e9) && (std::abs(b) <= 1e9);
    });
    valid = bad ? 0 : 1;
  }

  MPI_Bcast(&valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return valid == 1;
}

bool EreminVIntegralsMonteCarloALL::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool EreminVIntegralsMonteCarloALL::RunImpl() {
  const auto &input = GetInput();
  const auto &bounds = input.bounds;
  const int total_samples = input.samples;
  const auto &func = input.func;
  const std::size_t dimension = bounds.size();

  int rank = 0;
  int num_procs = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  double volume = 1.0;
  for (const auto &[a, b] : bounds) {
    volume *= (b - a);
  }

  const int base_samples = total_samples / num_procs;
  const int remainder = total_samples % num_procs;
  const int local_samples = base_samples + (rank == 0 ? remainder : 0);

  double local_sum = 0.0;

#pragma omp parallel reduction(+ : local_sum) default(none) \
    shared(bounds, local_samples, func, dimension, rank, num_procs)
  {
    const unsigned int seed = static_cast<unsigned int>(std::random_device{}()) +
                              (static_cast<unsigned int>(rank) * 1000U) +
                              static_cast<unsigned int>(omp_get_thread_num());

    std::mt19937 local_gen(seed);

    std::vector<std::uniform_real_distribution<double>> local_dists;
    local_dists.reserve(dimension);
    for (const auto &[a, b] : bounds) {
      local_dists.emplace_back(a, b);
    }

    std::vector<double> point(dimension);

#pragma omp for schedule(static)
    for (int i = 0; i < local_samples; ++i) {
      for (std::size_t dim = 0; dim < dimension; ++dim) {
        point[dim] = local_dists[dim](local_gen);
      }
      local_sum += func(point);
    }
  }

  double global_sum = 0.0;
  MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  double result = 0.0;
  if (rank == 0) {
    result = volume * (global_sum / static_cast<double>(total_samples));
  }
  MPI_Bcast(&result, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  GetOutput() = result;
  return true;
}

bool EreminVIntegralsMonteCarloALL::PostProcessingImpl() {
  return true;
}

}  // namespace eremin_v_integrals_monte_carlo
