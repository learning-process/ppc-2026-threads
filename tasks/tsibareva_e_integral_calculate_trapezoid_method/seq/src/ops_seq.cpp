#include "tsibareva_e_integral_calculate_trapezoid_method/seq/include/ops_seq.hpp"

#include <cmath>
#include <vector>

#include "task/include/task.hpp"
#include "tsibareva_e_integral_calculate_trapezoid_method/common/include/common.hpp"

namespace tsibareva_e_integral_calculate_trapezoid_method {

TsibarevaEIntegralCalculateTrapezoidMethodSEQ::TsibarevaEIntegralCalculateTrapezoidMethodSEQ(const Integral &in)
    : ppc::task::Task<Integral, double>() {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodSEQ::ValidationImpl() {
  return true;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodSEQ::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

std::vector<double> TsibarevaEIntegralCalculateTrapezoidMethodSEQ::ComputePoint(const std::vector<int> &indexes,
                                                                                const std::vector<double> &h, int dim) {
  std::vector<double> point(dim);
  for (int i = 0; i < dim; ++i) {
    point[i] = GetInput().lo[i] + (indexes[i] * h[i]);
  }
  return point;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodSEQ::IterateGridPoints(std::vector<int> &indexes, int dim) {
  int position = dim - 1;
  while (position >= 0) {
    indexes[position]++;
    if (indexes[position] <= GetInput().steps[position]) {
      return true;
    }
    indexes[position] = 0;
    position--;
  }
  return false;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodSEQ::RunImpl() {
  int dim = GetInput().dim;

  std::vector<double> h(dim);
  for (int i = 0; i < dim; ++i) {
    h[i] = (GetInput().hi[i] - GetInput().lo[i]) / GetInput().steps[i];
  }

  std::vector<int> indexes(dim, 0);
  double sum = 0.0;

  while (true) {
    std::vector<double> point = ComputePoint(indexes, h, dim);

    int boundary_count = 0;
    for (int i = 0; i < dim; ++i) {
      if (indexes[i] == 0 || indexes[i] == GetInput().steps[i]) {
        boundary_count++;
      }
    }

    double weight = (boundary_count == 0) ? 1.0 : std::pow(0.5, boundary_count);
    sum += weight * GetInput().f(point);

    if (!IterateGridPoints(indexes, dim)) {
      break;
    }
  }

  double res_h = 1.0;
  for (int i = 0; i < dim; ++i) {
    res_h *= h[i];
  }

  GetOutput() = sum * res_h;

  return true;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace tsibareva_e_integral_calculate_trapezoid_method
