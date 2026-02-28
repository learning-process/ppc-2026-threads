#include "samoylenko_i_integral_trapezoid/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

#include "samoylenko_i_integral_trapezoid/common/include/common.hpp"

namespace samoylenko_i_integral_trapezoid {

namespace {
std::function<double(const std::vector<double> &)> GetIntegrationFunction(int choice) {
  switch (choice) {
    case 0:  // x + y + ...
      return [](const std::vector<double> &values) {
        double sum = 0.0;
        for (double val : values) {
          sum += val;
        }
        return sum;
      };
    case 1:  // x * y * ...
      return [](const std::vector<double> &values) {
        double mult = 1.0;
        for (double val : values) {
          mult *= val;
        }
        return mult;
      };
    case 2:  // x^2 + y^2 + ...
      return [](const std::vector<double> &values) {
        double sum = 0.0;
        for (double val : values) {
          sum += val * val;
        }
        return sum;
      };
    case 3:  // sin(x + y + ...)
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

SamoylenkoITrapezoidIntegrationSEQ::SamoylenkoITrapezoidIntegrationSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SamoylenkoITrapezoidIntegrationSEQ::ValidationImpl() {
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

bool SamoylenkoITrapezoidIntegrationSEQ::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool SamoylenkoITrapezoidIntegrationSEQ::RunImpl() {
  const auto &in = GetInput();
  const int dimensions = static_cast<int>(in.a.size());
  auto integral_function = GetIntegrationFunction(in.function_choice);
  double sum = 0.0;

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

  std::vector<double> current_point(dimensions);

  for (int64_t pnt = 0; pnt < points; pnt++) {
    int64_t rem_index = pnt;
    int weight = 1;

    for (int dim = 0; dim < dimensions; dim++) {
      // convert point index to its coordinate in dimension "dim"
      int dim_coord = static_cast<int>(rem_index % dim_sizes[dim]);
      rem_index /= dim_sizes[dim];

      current_point[dim] = in.a[dim] + (dim_coord * h[dim]);

      if (dim_coord > 0 && dim_coord < in.n[dim]) {
        weight *= 2;
      }
    }

    sum += integral_function(current_point) * weight;
  }

  double h_mult = 1.0;
  for (int i = 0; i < dimensions; i++) {
    h_mult *= h[i];
  }

  GetOutput() = sum * (h_mult / std::pow(2.0, dimensions));  // h1 * h2 * ... * hn / 2^n

  return true;
}

bool SamoylenkoITrapezoidIntegrationSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace samoylenko_i_integral_trapezoid
