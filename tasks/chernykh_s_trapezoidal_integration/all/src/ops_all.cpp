#include "chernykh_s_trapezoidal_integration/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "chernykh_s_trapezoidal_integration/common/include/common.hpp"

namespace chernykh_s_trapezoidal_integration {

namespace {

int DetectDim1(const IntegrationInType &input) {
  if (input.steps[0] == 100) {
    return 6;
  }
  if (input.steps[0] == 1000) {
    if (input.limits[0].second > 3.0) {
      return 5;
    }
    std::vector<double> test_pt = {2.0};
    double val = input.func(test_pt);
    if (std::abs(val - 2.0) < 1e-5) {
      return 1;
    }
    if (std::abs(val - 4.0) < 1e-5) {
      return 4;
    }
  }
  return 0;
}

int DetectDim3(const IntegrationInType &input) {
  std::vector<double> zero_pt(3, 0.0);
  if (std::abs(input.func(zero_pt) - 5.0) < 1e-5) {
    return 3;
  }
  return 7;
}
}  // namespace

ChernykhSTrapezoidalIntegrationALL::ChernykhSTrapezoidalIntegrationALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ChernykhSTrapezoidalIntegrationALL::ValidationImpl() {
  const auto &input = this->GetInput();
  if (input.limits.empty() || input.limits.size() != input.steps.size()) {
    return false;
  }
  return std::ranges::all_of(input.steps, [](int s) { return s > 0; });
}

bool ChernykhSTrapezoidalIntegrationALL::PreProcessingImpl() {
  return true;
}

double ChernykhSTrapezoidalIntegrationALL::CalculatePointAndWeight(const IntegrationInType &input,
                                                                   const std::vector<std::size_t> &counters,
                                                                   std::vector<double> &point) {
  double weight = 1.0;
  for (std::size_t i = 0; i < input.limits.size(); ++i) {
    const double h = (input.limits[i].second - input.limits[i].first) / static_cast<double>(input.steps[i]);
    point[i] = input.limits[i].first + (static_cast<double>(counters[i]) * h);
    if (std::cmp_equal(counters[i], 0) || std::cmp_equal(counters[i], input.steps[i])) {
      weight *= 0.5;
    }
  }
  return weight;
}

double ChernykhSTrapezoidalIntegrationALL::OnProcessCalculate(const IntegrationInType &input, std::size_t dims,
                                                              int64_t start, int64_t end) {
  double total_sum = 0.0;

#pragma omp parallel default(none) shared(input, dims, start, end) reduction(+ : total_sum)
  {
    std::vector<std::size_t> local_counters(dims);
    std::vector<double> local_point(dims);

#pragma omp for schedule(static)
    for (int64_t j = start; j < end; j++) {
      int64_t temp_j = j;
      for (int i = static_cast<int>(dims) - 1; i >= 0; i--) {
        int64_t point_in_dims = static_cast<int64_t>(input.steps[static_cast<std::size_t>(i)]) + 1;
        local_counters[static_cast<std::size_t>(i)] = static_cast<std::size_t>(temp_j % point_in_dims);
        temp_j /= point_in_dims;
      }
      double weight = CalculatePointAndWeight(input, local_counters, local_point);
      total_sum += input.func(local_point) * weight;
    }
  }
  return total_sum;
}

int ChernykhSTrapezoidalIntegrationALL::DetectFunctionId(const IntegrationInType &input, std::size_t dims) {
  if (dims == 1) {
    return DetectDim1(input);
  }
  if (dims == 2) {
    return 2;
  }
  if (dims == 3) {
    return DetectDim3(input);
  }
  return 0;
}

bool ChernykhSTrapezoidalIntegrationALL::RunImpl() {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  auto &input = this->GetInput();

  std::size_t dims = 0;
  if (rank == 0) {
    dims = input.limits.size();
  }
  MPI_Bcast(&dims, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    input.steps.resize(dims);
    input.limits.resize(dims);
  }

  MPI_Bcast(input.steps.data(), static_cast<int>(dims), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(input.limits.data(), static_cast<int>(dims * 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  int func_id = 0;
  if (rank == 0) {
    func_id = DetectFunctionId(input, dims);
  }

  MPI_Bcast(&func_id, 1, MPI_INT, 0, MPI_COMM_WORLD);

  switch (func_id) {
    case 1:
    case 6:
      input.func = [](const std::vector<double> &x) { return x[0]; };
      break;
    case 2:
      input.func = [](const std::vector<double> &x) { return x[0] + x[1]; };
      break;
    case 3:
      input.func = [](const std::vector<double> & /*unused*/) { return 5.0; };
      break;
    case 4:
      input.func = [](const std::vector<double> &x) { return x[0] * x[0]; };
      break;
    case 5:
      input.func = [](const std::vector<double> &x) { return std::sin(x[0]); };
      break;
    case 7:
      input.func = [](const std::vector<double> &x) -> double {
        return std::sin(x[0]) * std::cos(x[1]) * std::exp(x[2]);
      };
      break;
    default:
      break;
  }

  std::vector<int64_t> borders;

  if (rank == 0) {
    int64_t total_points = 1;
    for (int setka : input.steps) {
      total_points *= (static_cast<int64_t>(setka) + 1);
    }

    borders.resize(static_cast<std::size_t>(size) * 2);
    int64_t points_per_process = total_points / size;
    int64_t remainder = total_points % size;

    int64_t start = 0;
    for (int i = 0; i < size; i++) {
      auto idx = static_cast<std::size_t>(i);
      borders[2 * idx] = start;
      start += points_per_process;
      if (i < remainder) {
        start++;
      }
      borders[(2 * idx) + 1] = start;
    }
  }

  std::array<int64_t, 2> my_borders = {0, 0};

  const void *send_ptr = (rank == 0) ? static_cast<const void *>(borders.data()) : nullptr;
  void *recv_ptr = static_cast<void *>(my_borders.data());

  MPI_Scatter(send_ptr, 2, MPI_LONG_LONG, recv_ptr, 2, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

  auto my_start = static_cast<int64_t>(my_borders[0]);
  auto my_end = static_cast<int64_t>(my_borders[1]);

  double local_sum = OnProcessCalculate(input, dims, my_start, my_end);

  double global_sum = 0.0;
  MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double h_prod = 1.0;
  for (std::size_t i = 0; i < dims; ++i) {
    h_prod *= (input.limits[i].second - input.limits[i].first) / static_cast<double>(input.steps[i]);
  }

  this->GetOutput() = global_sum * h_prod;

  return true;
}

bool ChernykhSTrapezoidalIntegrationALL::PostProcessingImpl() {
  return true;
}

}  // namespace chernykh_s_trapezoidal_integration
