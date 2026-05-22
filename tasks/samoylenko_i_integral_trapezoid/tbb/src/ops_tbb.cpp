#include "samoylenko_i_integral_trapezoid/tbb/include/ops_tbb.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_reduce.h"
#include "samoylenko_i_integral_trapezoid/common/include/common.hpp"

namespace samoylenko_i_integral_trapezoid {

namespace {
std::function<double(const std::vector<double> &)> GetIntegrationFunction(int64_t choice) {
  switch (choice) {
    case 0:
      return [](const std::vector<double> &values) {
        double sum = 0.0;
        for (double val : values) {
          sum += val;
        }
        return sum;
      };
    case 1:
      return [](const std::vector<double> &values) {
        double mult = 1.0;
        for (double val : values) {
          mult *= val;
        }
        return mult;
      };
    case 2:
      return [](const std::vector<double> &values) {
        double sum = 0.0;
        for (double val : values) {
          sum += val * val;
        }
        return sum;
      };
    case 3:
      return [](const std::vector<double> &values) {
        double sum = 0.0;
        for (double val : values) {
          sum += val;
        }
        return std::sin(sum);
      };
    default:
      return [](const std::vector<double> &) { return 0.0; };
  }
}
}  // namespace

SamoylenkoIIntegralTrapezoidTBB::SamoylenkoIIntegralTrapezoidTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SamoylenkoIIntegralTrapezoidTBB::ValidationImpl() {
  const auto &in = GetInput();

  if (in.a.empty() || in.a.size() != in.b.size() || in.a.size() != in.n.size()) {
    return false;
  }

  for (size_t i = 0; i < in.a.size(); ++i) {
    if (in.n[i] <= 0 || in.a[i] >= in.b[i]) {
      return false;
    }
  }

  return in.function_choice >= 0 && in.function_choice <= 3;
}

bool SamoylenkoIIntegralTrapezoidTBB::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool SamoylenkoIIntegralTrapezoidTBB::RunImpl() {
  const auto &in = GetInput();
  const int dimensions = static_cast<int>(in.a.size());
  auto integral_function = GetIntegrationFunction(in.function_choice);

  std::vector<double> h(dimensions);
  for (int i = 0; i < dimensions; i++) {
    h[i] = (in.b[i] - in.a[i]) / in.n[i];
  }

  std::vector<int64_t> dim_sizes(dimensions);
  int64_t points = 1;
  for (int i = 0; i < dimensions; i++) {
    dim_sizes[i] = in.n[i] + 1;
    points *= dim_sizes[i];
  }

  double sum = tbb::parallel_reduce(tbb::blocked_range<int64_t>(0, points), 0.0,
                                    [&h, &dimensions, &dim_sizes, &in, &integral_function](
                                        const tbb::blocked_range<int64_t> &range, double local_sum) -> double {
    std::vector<double> current_point(dimensions);

    for (int64_t pnt = range.begin(); pnt != range.end(); ++pnt) {
      int64_t rem_index = pnt;
      int weight = 1;

      for (int dim = 0; dim < dimensions; dim++) {
        int dim_coord = static_cast<int>(rem_index % dim_sizes[dim]);
        rem_index /= dim_sizes[dim];

        current_point[dim] = in.a[dim] + (dim_coord * h[dim]);

        if (dim_coord > 0 && dim_coord < in.n[dim]) {
          weight *= 2;
        }
      }

      local_sum += integral_function(current_point) * weight;
    }

    return local_sum;
  },
                                    [](double local_sum1, double local_sum2) { return local_sum1 + local_sum2; });

  double h_mult = 1.0;
  for (int i = 0; i < dimensions; i++) {
    h_mult *= h[i];
  }

  GetOutput() = sum * (h_mult / std::pow(2.0, dimensions));

  return true;
}

bool SamoylenkoIIntegralTrapezoidTBB::PostProcessingImpl() {
  return true;
}

}  // namespace samoylenko_i_integral_trapezoid
