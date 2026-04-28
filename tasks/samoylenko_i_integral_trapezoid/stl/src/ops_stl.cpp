#include "samoylenko_i_integral_trapezoid/stl/include/ops_stl.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <thread>
#include <vector>

#include "samoylenko_i_integral_trapezoid/common/include/common.hpp"
#include "util/include/util.hpp"

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

SamoylenkoIIntegralTrapezoidSTL::SamoylenkoIIntegralTrapezoidSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SamoylenkoIIntegralTrapezoidSTL::ValidationImpl() {
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

bool SamoylenkoIIntegralTrapezoidSTL::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

namespace {
double GetLocalSum(int64_t start, int64_t end, int dimensions, const std::vector<int64_t> &dim_sizes,
                   const std::vector<double> &h, const auto &in, auto &integral_function) {
  std::vector<double> current_point(dimensions);
  double local_sum = 0.0;

  for (int64_t pnt = start; pnt < end; pnt++) {
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
}
}  // namespace

bool SamoylenkoIIntegralTrapezoidSTL::RunImpl() {
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

  const int num_threads = ppc::util::GetNumThreads();
  std::vector<double> local_sums(num_threads, 0.0);
  std::vector<std::thread> threads(num_threads);

  for (int thr = 0; thr < num_threads; thr++) {
    int64_t start = (points * thr) / num_threads;
    int64_t end = (points * (thr + 1)) / num_threads;

    threads[thr] = std::thread([&h, &dimensions, &dim_sizes, &in, &integral_function, thr, start, end, &local_sums]() {
      local_sums[thr] = GetLocalSum(start, end, dimensions, dim_sizes, h, in, integral_function);
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }

  double sum = 0.0;
  for (int thr = 0; thr < num_threads; thr++) {
    sum += local_sums[thr];
  }

  double h_mult = 1.0;
  for (int i = 0; i < dimensions; i++) {
    h_mult *= h[i];
  }

  GetOutput() = sum * (h_mult / std::pow(2.0, dimensions));

  return true;
}

bool SamoylenkoIIntegralTrapezoidSTL::PostProcessingImpl() {
  return true;
}

}  // namespace samoylenko_i_integral_trapezoid
