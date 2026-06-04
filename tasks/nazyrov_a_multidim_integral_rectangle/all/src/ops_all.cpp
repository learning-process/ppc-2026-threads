#include "nazyrov_a_multidim_integral_rectangle/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "nazyrov_a_multidim_integral_rectangle/common/include/common.hpp"

namespace nazyrov_a_multidim_integral_rectangle {

NazyrovAMultidimIntegralRectangleAll::NazyrovAMultidimIntegralRectangleAll(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool NazyrovAMultidimIntegralRectangleAll::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) {
    return true;
  }
  const auto &func = std::get<0>(GetInput());
  const auto &bounds = std::get<1>(GetInput());
  const int n = std::get<2>(GetInput());
  return func && n > 0 && !bounds.empty() && std::ranges::all_of(bounds, [](const auto &bd) {
    return std::isfinite(bd.first) && std::isfinite(bd.second) && bd.first < bd.second;
  });
}

bool NazyrovAMultidimIntegralRectangleAll::PreProcessingImpl() {
  return true;
}

bool NazyrovAMultidimIntegralRectangleAll::RunImpl() {
  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  const auto &func = std::get<0>(GetInput());
  const auto &bounds = std::get<1>(GetInput());
  const int n = std::get<2>(GetInput());

  const int dim = static_cast<int>(bounds.size());

  std::vector<double> h(static_cast<std::size_t>(dim));
  double cell_vol = 1.0;
  for (int i = 0; i < dim; ++i) {
    h[i] = (bounds[i].second - bounds[i].first) / n;
    cell_vol *= h[i];
  }

  std::int64_t total = 1;
  for (int i = 0; i < dim; ++i) {
    total *= n;
  }

  const auto rk = static_cast<std::int64_t>(rank);
  const auto ws = static_cast<std::int64_t>(world_size);
  const std::int64_t cells_per_rank = total / ws;
  const std::int64_t extra = total % ws;
  const std::int64_t begin = (rk * cells_per_rank) + std::min(rk, extra);
  const std::int64_t end = begin + cells_per_rank + (rk < extra ? 1LL : 0LL);

  double local_sum = 0.0;

#pragma omp parallel default(none) shared(func, bounds, h, n, dim, begin, end, local_sum)
  {
    std::vector<double> point(static_cast<std::size_t>(dim));
    double thread_sum = 0.0;

#pragma omp for schedule(static)
    for (std::int64_t cell = begin; cell < end; ++cell) {
      std::int64_t tmp = cell;
      for (int i = dim - 1; i >= 0; --i) {
        const int ki = static_cast<int>(tmp % n);
        tmp /= n;
        const double coordinate = bounds[i].first + ((ki + 0.5) * h[i]);
        point[i] = coordinate;
      }
      thread_sum += func(point);
    }

#pragma omp atomic
    local_sum += thread_sum;
  }

  double global_sum = 0.0;
  MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    GetOutput() = global_sum * cell_vol;
  }

  return true;
}

bool NazyrovAMultidimIntegralRectangleAll::PostProcessingImpl() {
  return true;
}

}  // namespace nazyrov_a_multidim_integral_rectangle
