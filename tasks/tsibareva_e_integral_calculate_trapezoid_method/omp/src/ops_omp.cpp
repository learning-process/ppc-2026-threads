#include "tsibareva_e_integral_calculate_trapezoid_method/omp/include/ops_omp.hpp"

#include <cmath>
#include <vector>

#include "tsibareva_e_integral_calculate_trapezoid_method/common/include/common.hpp"

namespace tsibareva_e_integral_calculate_trapezoid_method {

TsibarevaEIntegralCalculateTrapezoidMethodOMP::TsibarevaEIntegralCalculateTrapezoidMethodOMP(const InType &in)
    : BaseTask() {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodOMP::ValidationImpl() {
  return true;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodOMP::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodOMP::RunImpl() {
  int dim = GetInput().dim;

  std::vector<double> h(dim);
  std::vector<int> sizes(dim);
  int total_nodes = 1;
  for (int i = 0; i < dim; ++i) {
    h[i] = (GetInput().hi[i] - GetInput().lo[i]) / GetInput().steps[i];
    sizes[i] = GetInput().steps[i] + 1;
    total_nodes *= sizes[i];
  }

  double global_sum = 0.0;

#pragma omp parallel for default(none) shared(dim, h, sizes, total_nodes) reduction(+ : global_sum)
  for (int node = 0; node < total_nodes; ++node) {
    int remainder = node;
    double node_weight = 1.0;
    std::vector<double> point(dim);

    for (int i = dim - 1; i >= 0; --i) {
      int idx = remainder % sizes[i];
      remainder /= sizes[i];

      if (idx == 0 || idx == GetInput().steps[i]) {
        node_weight *= 0.5;
      }

      point[i] = GetInput().lo[i] + (idx * h[i]);
    }

    global_sum += node_weight * GetInput().f(point);
  }

  double res_h = 1.0;
  for (int i = 0; i < dim; ++i) {
    res_h *= h[i];
  }
  GetOutput() = global_sum * res_h;
  return true;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodOMP::PostProcessingImpl() {
  return true;
}

}  // namespace tsibareva_e_integral_calculate_trapezoid_method
