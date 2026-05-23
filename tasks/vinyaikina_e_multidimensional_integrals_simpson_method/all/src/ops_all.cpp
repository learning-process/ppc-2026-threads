#include "vinyaikina_e_multidimensional_integrals_simpson_method/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

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

VinyaikinaEMultidimIntegrSimpsonALL::VinyaikinaEMultidimIntegrSimpsonALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool VinyaikinaEMultidimIntegrSimpsonALL::PreProcessingImpl() {
  return true;
}

bool VinyaikinaEMultidimIntegrSimpsonALL::ValidationImpl() {
  const auto &[h, limits, function] = GetInput();
  return !limits.empty() && function && h <= 0.01;
}

bool VinyaikinaEMultidimIntegrSimpsonALL::RunImpl() {
  const auto &input = GetInput();
  double h = std::get<0>(input);
  const auto &limits = std::get<1>(input);
  const auto &function = std::get<2>(input);

  const int num_threads = ppc::util::GetNumThreads();

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

  int mpi_rank = 0;
  int mpi_size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  double delta_mpi = (limits[0].second - limits[0].first) / mpi_size;

  double proc_left = limits[0].first;
  double proc_right = limits[0].second;

  if (mpi_rank != 0) {
    proc_left = CustomRound(limits[0].first + (delta_mpi * mpi_rank), actual_step[0]);
  }
  if (mpi_rank != mpi_size - 1) {
    proc_right = CustomRound(limits[0].second - (delta_mpi * (mpi_size - mpi_rank - 1)), actual_step[0]);
  }

  double res = 0.0;
  double delta_omp = (proc_right - proc_left) / num_threads;

#pragma omp parallel num_threads(num_threads) default(none)                                              \
    shared(proc_left, proc_right, limits, simpson_factor, actual_step, delta_omp, function, num_threads) \
    reduction(+ : res)
  {
    double left_border = proc_left;
    double right_border = proc_right;
    int tid = omp_get_thread_num();

    if (tid != 0) {
      left_border = CustomRound(proc_left + (delta_omp * tid), actual_step[0]);
    }
    if (tid != num_threads - 1) {
      right_border = CustomRound(proc_right - (delta_omp * (num_threads - tid - 1)), actual_step[0]);
    }

    res += OuntNtIntegral(left_border, right_border, simpson_factor, limits, actual_step, function);
  }

  double global_res = 0.0;
  MPI_Reduce(&res, &global_res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(&global_res, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  I_res_ = global_res;

  return true;
}

bool VinyaikinaEMultidimIntegrSimpsonALL::PostProcessingImpl() {
  GetOutput() = I_res_;
  return true;
}
}  // namespace vinyaikina_e_multidimensional_integrals_simpson_method
