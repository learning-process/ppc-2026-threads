#include "galkin_d_multidim_integrals_rectangles/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <tuple>
#include <vector>

#include "galkin_d_multidim_integrals_rectangles/common/include/common.hpp"
#include "util/include/util.hpp"

namespace galkin_d_multidim_integrals_rectangles {

GalkinDMultidimIntegralsRectanglesALL::GalkinDMultidimIntegralsRectanglesALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool GalkinDMultidimIntegralsRectanglesALL::ValidationImpl() {
  const auto &func = std::get<0>(GetInput());
  const auto &borders = std::get<1>(GetInput());
  const int n = std::get<2>(GetInput());

  if (borders.empty()) {
    return false;
  }
  for (const auto &[left_border, right_border] : borders) {
    if (!std::isfinite(left_border) || !std::isfinite(right_border)) {
      return false;
    }
    if (left_border >= right_border) {
      return false;
    }
  }
  return func && (n > 0) && (GetOutput() == 0.0);
}

bool GalkinDMultidimIntegralsRectanglesALL::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool GalkinDMultidimIntegralsRectanglesALL::RunImpl() {
  const InType &input = GetInput();
  const auto &func = std::get<0>(input);
  const auto &borders = std::get<1>(input);
  const int n = std::get<2>(input);
  const std::size_t dim = borders.size();

  std::vector<double> h(dim);
  double cell_v = 1.0;
  for (std::size_t i = 0; i < dim; ++i) {
    h[i] = (borders[i].second - borders[i].first) / static_cast<double>(n);
    if (!(h[i] > 0.0) || !std::isfinite(h[i])) {
      return false;
    }
    cell_v *= h[i];
  }

  std::size_t total = 1;
  for (std::size_t i = 0; i < dim; ++i) {
    if (total > (std::numeric_limits<std::size_t>::max() / static_cast<std::size_t>(n))) {
      return false;
    }
    total *= static_cast<std::size_t>(n);
  }
  if (total > static_cast<std::size_t>(LLONG_MAX)) {
    return false;
  }
  const auto total_i64 = static_cast<std::int64_t>(total);

  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  const std::int64_t chunk = total_i64 / world_size;
  const std::int64_t rem = total_i64 % world_size;
  const std::int64_t mpi_begin =
      (static_cast<std::int64_t>(rank) * chunk) + std::min(static_cast<std::int64_t>(rank), rem);
  const std::int64_t mpi_end = mpi_begin + chunk + (static_cast<std::int64_t>(rank) < rem ? 1 : 0);

  double local_sum = 0.0;

#pragma omp parallel for reduction(+ : local_sum) schedule(static)     \
    num_threads(std::max(1, ppc::util::GetNumThreads())) default(none) \
    shared(borders, h, dim, func, n, mpi_begin, mpi_end)
  for (std::int64_t linear_idx = mpi_begin; linear_idx < mpi_end; ++linear_idx) {
    std::vector<double> x(dim);
    auto tmp = static_cast<std::size_t>(linear_idx);
    for (std::size_t i = 0; i < dim; ++i) {
      const std::size_t idx_i = tmp % static_cast<std::size_t>(n);
      tmp /= static_cast<std::size_t>(n);
      x[i] = borders[i].first + ((static_cast<double>(idx_i) + 0.5) * h[i]);
    }
    local_sum += func(x);
  }

  double global_sum = 0.0;
  MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(&global_sum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  GetOutput() = global_sum * cell_v;
  return std::isfinite(GetOutput());
}

bool GalkinDMultidimIntegralsRectanglesALL::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace galkin_d_multidim_integrals_rectangles
