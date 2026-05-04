#include "vinyaikina_e_multidimensional_integrals_simpson_method/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <stack>
#include <utility>
#include <vector>

#include "vinyaikina_e_multidimensional_integrals_simpson_method/common/include/common.hpp"

namespace vinyaikina_e_multidimensional_integrals_simpson_method {
namespace {

double OuntNtIntegral(double simpson_factor, const std::vector<std::pair<double, double>> &limits,
                      std::vector<double> &actual_step,
                      const std::function<double(const std::vector<double> &)> &function) {
  std::stack<std::pair<std::vector<double>, double>> stack;
  stack.emplace(std::vector<double>(), 1.0);

  double res = 0.0;

  while (!stack.empty()) {
    std::vector<double> point = stack.top().first;
    double weight = stack.top().second;
    stack.pop();

    if (point.size() == limits.size()) {
      res += function(point) * weight * simpson_factor;
      continue;
    }

    size_t dim = point.size();
    double step = actual_step[dim] / 1.0;

    int steps_count = static_cast<int>(lround((limits[dim].second - limits[dim].first) / step));

    for (int i = 0; i <= steps_count; ++i) {
      double x = limits[dim].first + (i * step);

      double dim_weight = 2.0;
      if (i == 0 || i == steps_count) {
        dim_weight = 1.0;
      } else if (i % 2 != 0) {
        dim_weight = 4.0;
      }

      point.push_back(x);
      stack.emplace(point, weight * dim_weight);
      point.pop_back();
    }
  }

  return res;
}
};  // namespace

VinyaikinaEMultidimIntegrSimpsonSEQ::VinyaikinaEMultidimIntegrSimpsonSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool VinyaikinaEMultidimIntegrSimpsonSEQ::PreProcessingImpl() {
  return true;
}

bool VinyaikinaEMultidimIntegrSimpsonSEQ::ValidationImpl() {
  const auto &[h, limits, function] = GetInput();
  return !limits.empty() && function && h <= 0.01;
}

bool VinyaikinaEMultidimIntegrSimpsonSEQ::RunImpl() {
  const auto &[h, limits, function] = GetInput();

  std::vector<double> actual_step(limits.size());
  double simpson_factor = 1.0;

  for (size_t i = 0; i < limits.size(); i++) {
    int quan_steps = static_cast<int>(lround((limits[i].second - limits[i].first) / h));
    if (quan_steps % 2 != 0) {
      quan_steps++;
    }
    actual_step[i] = (limits[i].second - limits[i].first) / quan_steps;
    simpson_factor *= actual_step[i] / 3.0;
  }

  I_res_ = OuntNtIntegral(simpson_factor, limits, actual_step, function);

  return true;
}

bool VinyaikinaEMultidimIntegrSimpsonSEQ::PostProcessingImpl() {
  GetOutput() = I_res_;
  return true;
}
}  // namespace vinyaikina_e_multidimensional_integrals_simpson_method
