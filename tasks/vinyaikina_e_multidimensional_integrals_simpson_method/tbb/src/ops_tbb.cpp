#include "vinyaikina_e_multidimensional_integrals_simpson_method/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <stack>
#include <utility>
#include <vector>

#include "util/include/util.hpp"
#include "vinyaikina_e_multidimensional_integrals_simpson_method/common/include/common.hpp"

namespace vinyaikina_e_multidimensional_integrals_simpson_method {
namespace {

double CustomRound(double value, double h) {
  h *= 2;
  int tmp = static_cast<int>(1 / h);
  int decimal_places = 0;
  while (tmp > 0 && tmp % 10 == 0) {
    decimal_places++;
    tmp /= 10;
  }

  double factor = std::pow(10.0, decimal_places);
  return std::round(value * factor) / factor;
}

double Weight(int i, int steps_count) {
  double weight = 2.0;
  if (i == 0 || i == steps_count) {
    weight = 1.0;
  } else if (i % 2 != 0) {
    weight = 4.0;
  }
  return weight;
}

double OuntNtIntegral(double left_border, double right_border, double simpson_factor,
                      const std::vector<std::pair<double, double>> &limits, const std::vector<double> &actual_step,
                      const std::function<double(const std::vector<double> &)> &function) {
  std::stack<std::pair<std::vector<double>, double>> stack;
  double res = 0.0;

  int steps_count_0 = static_cast<int>(lround((right_border - left_border) / actual_step[0]));

  for (int i0 = 0; i0 <= steps_count_0; ++i0) {
    double x0 = left_border + (i0 * actual_step[0]);

    double weight_0 = Weight(i0, steps_count_0);
    stack.emplace(std::vector<double>{x0}, weight_0);

    while (!stack.empty()) {
      std::vector<double> point = stack.top().first;
      double weight = stack.top().second;
      stack.pop();

      if (point.size() == limits.size()) {
        res += function(point) * weight * simpson_factor;
        continue;
      }

      size_t dim = point.size();
      double step = actual_step[dim];

      int steps_count = static_cast<int>(lround((limits[dim].second - limits[dim].first) / step));

      for (int i = 0; i <= steps_count; ++i) {
        double x = limits[dim].first + (i * step);

        double dim_weight = Weight(i, steps_count);

        point.push_back(x);
        stack.emplace(point, weight * dim_weight);
        point.pop_back();
      }
    }
  }

  return res;
}
};  // namespace

VinyaikinaEMultidimIntegrSimpsonTBB::VinyaikinaEMultidimIntegrSimpsonTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool VinyaikinaEMultidimIntegrSimpsonTBB::PreProcessingImpl() {
  return true;
}

bool VinyaikinaEMultidimIntegrSimpsonTBB::ValidationImpl() {
  const auto &[h, limits, function] = GetInput();
  return !limits.empty() && function && h <= 0.01;
}

bool VinyaikinaEMultidimIntegrSimpsonTBB::RunImpl() {
  const auto &input = GetInput();
  double h = std::get<0>(input);
  const auto &limits = std::get<1>(input);
  const auto &function = std::get<2>(input);

  const int num_threads = ppc::util::GetNumThreads();

  double delta = (limits[0].second - limits[0].first) / num_threads;

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

  tbb::blocked_range<int> range(0, num_threads);

  I_res_ = tbb::parallel_reduce(range, 0.0, [&](const tbb::blocked_range<int> &r, double local_res) {
    for (int i = r.begin(); i != r.end(); ++i) {
      double left_border = limits[0].first;
      double right_border = limits[0].second;

      if (i != 0) {
        left_border = CustomRound(limits[0].first + (delta * i), actual_step[0]);
      }
      if (i != num_threads - 1) {
        right_border = CustomRound(limits[0].second - (delta * (num_threads - i - 1)), actual_step[0]);
      }

      local_res += OuntNtIntegral(left_border, right_border, simpson_factor, limits, actual_step, function);
    }
    return local_res;
  }, std::plus<>());

  return true;
}

bool VinyaikinaEMultidimIntegrSimpsonTBB::PostProcessingImpl() {
  GetOutput() = I_res_;
  return true;
}
}  // namespace vinyaikina_e_multidimensional_integrals_simpson_method
