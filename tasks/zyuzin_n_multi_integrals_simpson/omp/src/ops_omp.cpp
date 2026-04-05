#include "zyuzin_n_multi_integrals_simpson/omp/include/ops_omp.hpp"

#include <cstddef>
#include <vector>

#include "zyuzin_n_multi_integrals_simpson/common/include/common.hpp"

namespace zyuzin_n_multi_integrals_simpson {

ZyuzinNSimpsonOMP::ZyuzinNSimpsonOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ZyuzinNSimpsonOMP::ValidationImpl() {
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

bool ZyuzinNSimpsonOMP::PreProcessingImpl() {
  result_ = 0.0;
  return true;
}

double ZyuzinNSimpsonOMP::GetSimpsonWeight(int index, int n) {
  if (index == 0 || index == n) {
    return 1.0;
  }
  if (index % 2 == 1) {
    return 4.0;
  }
  return 2.0;
}

double ZyuzinNSimpsonOMP::ComputeSimpsonMultiDim() {
  const auto &input = GetInput();
  const size_t num_dims = input.lower_bounds.size();

  std::vector<double> h(num_dims);
  for (size_t dim = 0; dim < num_dims; ++dim) {
    h[dim] = (input.upper_bounds[dim] - input.lower_bounds[dim]) / input.n_steps[dim];
  }

  size_t total_points = 1;
  for (size_t dim = 0; dim < num_dims; ++dim) {
    total_points *= static_cast<size_t>(input.n_steps[dim] + 1);
  }

  double sum = 0.0;

#pragma omp parallel for default(none) shared(input, total_points, num_dims, h) reduction(+ : sum)
  for (size_t point_idx = 0; point_idx < total_points; ++point_idx) {
    std::vector<int> indices(num_dims);
    auto temp = point_idx;
    for (size_t dim = 0; dim < num_dims; ++dim) {
      indices[dim] = static_cast<int>(temp % static_cast<size_t>(input.n_steps[dim] + 1));
      temp /= static_cast<size_t>(input.n_steps[dim] + 1);
    }

    std::vector<double> point(num_dims);
    for (size_t dim = 0; dim < num_dims; ++dim) {
      point[dim] = input.lower_bounds[dim] + (indices[dim] * h[dim]);
    }

    double weight = 1.0;
    for (size_t dim = 0; dim < num_dims; ++dim) {
      weight *= GetSimpsonWeight(indices[dim], input.n_steps[dim]);
    }

    sum += weight * input.func(point);
  }

  double factor = 1.0;
  for (size_t dim = 0; dim < num_dims; ++dim) {
    factor *= h[dim] / 3.0;
  }

  return sum * factor;
}

bool ZyuzinNSimpsonOMP::RunImpl() {
  result_ = ComputeSimpsonMultiDim();
  return true;
}

bool ZyuzinNSimpsonOMP::PostProcessingImpl() {
  GetOutput() = result_;
  return true;
}

}  // namespace zyuzin_n_multi_integrals_simpson
