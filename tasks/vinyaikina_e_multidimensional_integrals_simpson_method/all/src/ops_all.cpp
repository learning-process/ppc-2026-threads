#include "vinyaikina_e_multidimensional_integrals_simpson_method/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
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

}  // namespace

VinyaikinaEMultidimIntegrSimpsonMPI::VinyaikinaEMultidimIntegrSimpsonMPI(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool VinyaikinaEMultidimIntegrSimpsonMPI::PreProcessingImpl() {
  return true;
}

bool VinyaikinaEMultidimIntegrSimpsonMPI::ValidationImpl() {
  const auto &[h, limits, function] = GetInput();
  return !limits.empty() && function && h <= 0.01;
}

bool VinyaikinaEMultidimIntegrSimpsonMPI::RunImpl() {
  const auto &input = GetInput();
  double h = std::get<0>(input);
  const auto &limits = std::get<1>(input);
  const auto &function = std::get<2>(input);

  int rank = 0, num_procs = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  int num_omp_threads = ppc::util::GetNumThreads();

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

  int total_steps_0 = static_cast<int>(lround((limits[0].second - limits[0].first) / actual_step[0]));
  int base_steps = total_steps_0 / num_procs;
  int remainder = total_steps_0 % num_procs;

  int my_start_step = rank * base_steps + std::min(rank, remainder);
  int my_end_step = (rank + 1) * base_steps + std::min(rank + 1, remainder);

  double mpi_left = limits[0].first + my_start_step * actual_step[0];
  double mpi_right = limits[0].first + my_end_step * actual_step[0];

  double omp_delta = (mpi_right - mpi_left) / num_omp_threads;
  double local_res = 0.0;

#pragma omp parallel num_threads(num_omp_threads) default(none)                                            \
    shared(mpi_left, mpi_right, omp_delta, num_omp_threads, actual_step, simpson_factor, limits, function) \
    reduction(+ : local_res)
  {
    int tid = omp_get_thread_num();
    double t_left = mpi_left + tid * omp_delta;
    double t_right = mpi_right - (num_omp_threads - tid - 1) * omp_delta;

    if (tid == 0) {
      t_left = mpi_left;
    }
    if (tid == num_omp_threads - 1) {
      t_right = mpi_right;
    }

    t_left = CustomRound(t_left, actual_step[0]);
    t_right = CustomRound(t_right, actual_step[0]);

    if (t_right > t_left) {
      local_res += OuntNtIntegral(t_left, t_right, simpson_factor, limits, actual_step, function);
    }
  }

  double global_res = 0.0;
  MPI_Reduce(&local_res, &global_res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    I_res_ = global_res;
  }
  return true;
}

bool VinyaikinaEMultidimIntegrSimpsonMPI::PostProcessingImpl() {
  GetOutput() = I_res_;
  return true;
}

}  // namespace vinyaikina_e_multidimensional_integrals_simpson_method
