#include "zyuzin_n_multi_integrals_simpson/all/include/ops_all.hpp"

#include <mpi.h>

#include <cstddef>
#include <vector>

#include "zyuzin_n_multi_integrals_simpson/common/include/common.hpp"

namespace zyuzin_n_multi_integrals_simpson {

ZyuzinNSimpsonALL::ZyuzinNSimpsonALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ZyuzinNSimpsonALL::ValidationImpl() {
  const auto &input = GetInput();
  if (input.lower_bounds.size() != input.upper_bounds.size() || input.lower_bounds.size() != input.n_steps.size()) {
    return false;
  }
  if (input.lower_bounds.empty()) {
    return false;
  }
  for (size_t i = 0; i < input.lower_bounds.size(); ++i) {
    if (input.lower_bounds[i] > input.upper_bounds[i]) {
      return false;
    }
    if (input.n_steps[i] <= 0 || input.n_steps[i] % 2 != 0) {
      return false;
    }
  }
  return static_cast<bool>(input.func);
}

bool ZyuzinNSimpsonALL::PreProcessingImpl() {
  result_ = 0.0;
  return true;
}

double ZyuzinNSimpsonALL::GetSimpsonWeight(int index, int n) {
  if (index == 0 || index == n) {
    return 1.0;
  }
  if (index % 2 == 1) {
    return 4.0;
  }
  return 2.0;
}

double ZyuzinNSimpsonALL::ComputeLocalSimpsonSum(size_t begin, size_t end) {
  const auto &input = GetInput();
  const size_t num_dims = input.lower_bounds.size();

  std::vector<double> h(num_dims);
  for (size_t dim = 0; dim < num_dims; ++dim) {
    h[dim] = (input.upper_bounds[dim] - input.lower_bounds[dim]) / static_cast<double>(input.n_steps[dim]);
  }

  double local_sum = 0.0;

#pragma omp parallel for default(none) shared(input, h, begin, end, num_dims) reduction(+ : local_sum)
  for (size_t point_idx = begin; point_idx < end; ++point_idx) {
    auto temp = point_idx;
    std::vector<double> point(num_dims);
    double weight = 1.0;

    for (size_t dim = 0; dim < num_dims; ++dim) {
      const auto axis_points = static_cast<size_t>(input.n_steps[dim]) + 1U;
      const auto index = static_cast<int>(temp % axis_points);
      temp /= axis_points;
      point[dim] = input.lower_bounds[dim] + (static_cast<double>(index) * h[dim]);
      weight *= GetSimpsonWeight(index, input.n_steps[dim]);
    }

    local_sum += weight * input.func(point);
  }

  return local_sum;
}

bool ZyuzinNSimpsonALL::RunImpl() {
  const auto &input = GetInput();
  const size_t num_dims = input.lower_bounds.size();

  size_t total_points = 1;
  for (size_t dim = 0; dim < num_dims; ++dim) {
    total_points *= static_cast<size_t>(input.n_steps[dim]) + 1U;
  }

  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const auto usize = static_cast<size_t>(size);
  const auto urank = static_cast<size_t>(rank);
  const size_t begin = (total_points * urank) / usize;
  const size_t end = (total_points * (urank + 1U)) / usize;

  const double local_sum = ComputeLocalSimpsonSum(begin, end);

  double global_sum = 0.0;
  MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  std::vector<double> h(num_dims);
  for (size_t dim = 0; dim < num_dims; ++dim) {
    h[dim] = (input.upper_bounds[dim] - input.lower_bounds[dim]) / static_cast<double>(input.n_steps[dim]);
  }

  double factor = 1.0;
  for (size_t dim = 0; dim < num_dims; ++dim) {
    factor *= h[dim] / 3.0;
  }

  result_ = global_sum * factor;
  return true;
}

bool ZyuzinNSimpsonALL::PostProcessingImpl() {
  GetOutput() = result_;
  return true;
}

}  // namespace zyuzin_n_multi_integrals_simpson
