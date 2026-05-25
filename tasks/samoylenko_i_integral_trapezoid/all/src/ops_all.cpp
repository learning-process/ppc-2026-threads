#include "samoylenko_i_integral_trapezoid/all/include/ops_all.hpp"

#include <mpi.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
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

double GetLocalSum(int64_t start, int64_t end, int dimensions, const std::vector<int64_t> &dim_sizes,
                   const std::vector<double> &h, const auto &in, auto &integral_function) {
  std::vector<double> current_point(dimensions);
  std::vector<int> coords(dimensions);
  double local_sum = 0.0;

  int64_t rem_index = start;
  for (int dim = 0; dim < dimensions; dim++) {
    coords[dim] = static_cast<int>(rem_index % dim_sizes[dim]);
    rem_index /= dim_sizes[dim];
  }

  for (int64_t pnt = start; pnt < end; pnt++) {
    int weight = 1;

    for (int dim = 0; dim < dimensions; dim++) {
      current_point[dim] = in.a[dim] + (coords[dim] * h[dim]);

      if (coords[dim] > 0 && coords[dim] < in.n[dim]) {
        weight *= 2;
      }
    }

    local_sum += integral_function(current_point) * weight;

    for (int dim = 0; dim < dimensions; dim++) {
      coords[dim]++;
      if (coords[dim] < dim_sizes[dim]) {
        break;
      }
      coords[dim] = 0;
    }
  }

  return local_sum;
}
}  // namespace

SamoylenkoIIntegralTrapezoidALL::SamoylenkoIIntegralTrapezoidALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SamoylenkoIIntegralTrapezoidALL::ValidationImpl() {
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

bool SamoylenkoIIntegralTrapezoidALL::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool SamoylenkoIIntegralTrapezoidALL::RunImpl() {
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

  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int64_t start = (points * rank) / size;
  int64_t end = (points * (rank + 1)) / size;

  const int num_threads = ppc::util::GetNumThreads();
  std::vector<double> local_sums(num_threads, 0.0);

#pragma omp parallel for num_threads(num_threads) default(none) \
    shared(num_threads, start, end, dimensions, dim_sizes, h, in, integral_function, local_sums)
  for (int thr = 0; thr < num_threads; thr++) {
    int64_t local_start = start + (((end - start) * thr) / num_threads);
    int64_t local_end = start + (((end - start) * (thr + 1)) / num_threads);
    local_sums[thr] = GetLocalSum(local_start, local_end, dimensions, dim_sizes, h, in, integral_function);
  }

  double sum = 0.0;
  for (int thr = 0; thr < num_threads; thr++) {
    sum += local_sums[thr];
  }

  double total_sum = 0.0;
  MPI_Allreduce(&sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double h_mult = 1.0;
  for (int i = 0; i < dimensions; i++) {
    h_mult *= h[i];
  }

  GetOutput() = total_sum * (h_mult / std::pow(2.0, dimensions));

  return true;
}

bool SamoylenkoIIntegralTrapezoidALL::PostProcessingImpl() {
  return true;
}

}  // namespace samoylenko_i_integral_trapezoid
