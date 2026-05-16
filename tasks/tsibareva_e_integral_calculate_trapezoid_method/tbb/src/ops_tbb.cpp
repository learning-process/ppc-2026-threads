// ops_tbb.cpp
#include "tsibareva_e_integral_calculate_trapezoid_method/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

#include <cmath>
#include <functional>
#include <vector>

#include "tsibareva_e_integral_calculate_trapezoid_method/common/include/common.hpp"

namespace tsibareva_e_integral_calculate_trapezoid_method {

TsibarevaEIntegralCalculateTrapezoidMethodTBB::TsibarevaEIntegralCalculateTrapezoidMethodTBB(const InType &in)
    : BaseTask() {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodTBB::ValidationImpl() {
  return true;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodTBB::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodTBB::RunImpl() {
  const Integral &input_i = GetInput();
  int dim = input_i.dim;

  std::vector<double> h(dim);
  std::vector<int> sizes(dim);
  int total_nodes = 1;
  for (int i = 0; i < dim; ++i) {
    h[i] = (input_i.hi[i] - input_i.lo[i]) / input_i.steps[i];
    sizes[i] = input_i.steps[i] + 1;
    total_nodes *= sizes[i];
  }

  double res_sum = tbb::parallel_reduce(tbb::blocked_range<int>(0, total_nodes), 0.0,
                                        [&](const tbb::blocked_range<int> &r, double local_sum) {
    return local_sum + ComputeSumNode(r, 0.0, input_i, h, sizes);
  }, std::plus<>());

  double res_h = 1.0;
  for (int i = 0; i < dim; ++i) {
    res_h *= h[i];
  }
  GetOutput() = res_sum * res_h;
  return true;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodTBB::PostProcessingImpl() {
  return true;
}

double TsibarevaEIntegralCalculateTrapezoidMethodTBB::ComputeSumNode(const tbb::blocked_range<int> &range, double init,
                                                                     const Integral &input_i,
                                                                     const std::vector<double> &h,
                                                                     const std::vector<int> &sizes) {
  double sum = init;
  int dim = input_i.dim;

  for (int node = range.begin(); node != range.end(); ++node) {
    int rem = node;
    double node_w = 1.0;
    std::vector<double> point(dim);

    for (int i = dim - 1; i >= 0; --i) {
      int idx = rem % sizes[i];
      rem /= sizes[i];

      if (idx == 0 || idx == input_i.steps[i]) {
        node_w *= 0.5;
      }

      point[i] = input_i.lo[i] + (idx * h[i]);
    }

    sum += node_w * input_i.f(point);
  }

  return sum;
}

}  // namespace tsibareva_e_integral_calculate_trapezoid_method
