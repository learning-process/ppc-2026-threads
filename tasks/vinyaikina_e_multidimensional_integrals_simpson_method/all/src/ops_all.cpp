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

double Weight(int i, int steps_count) {
  double weight = 2.0;
  if (i == 0 || i == steps_count) {
    weight = 1.0;
  } else if (i % 2 != 0) {
    weight = 4.0;
  }
  return weight;
}

double IntegrateRemainingDims(double x0, double simpson_factor, const std::vector<std::pair<double, double>> &limits,
                              const std::vector<double> &actual_step,
                              const std::function<double(const std::vector<double> &)> &function) {
  std::stack<std::pair<std::vector<double>, double>> stack;
  double result = 0.0;

  stack.emplace(std::vector<double>{x0}, 1.0);

  while (!stack.empty()) {
    auto [point, weight_product] = stack.top();
    stack.pop();

    size_t dim = point.size();

    if (dim == limits.size()) {
      result += function(point) * weight_product * simpson_factor;
      continue;
    }

    double step = actual_step[dim];
    int steps_count = static_cast<int>(lround((limits[dim].second - limits[dim].first) / step));
    for (int i = 0; i <= steps_count; ++i) {
      double x = limits[dim].first + i * step;
      double dim_weight = Weight(i, steps_count);
      point.push_back(x);
      stack.emplace(point, weight_product * dim_weight);
      point.pop_back();
    }
  }

  return result;
}

}  // namespace

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

  int rank = 0, num_procs = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  int num_omp_threads = ppc::util::GetNumThreads();

  std::vector<double> actual_step(limits.size());
  double simpson_factor = 1.0;
  for (size_t i = 0; i < limits.size(); ++i) {
    double len = limits[i].second - limits[i].first;
    int quan_steps = static_cast<int>(lround(len / h));
    if (quan_steps % 2 != 0) {
      ++quan_steps;
    }
    actual_step[i] = len / quan_steps;
    simpson_factor *= actual_step[i] / 3.0;
  }

  int total_steps_0 = static_cast<int>(lround((limits[0].second - limits[0].first) / actual_step[0]));

  int steps_per_proc = total_steps_0 / num_procs;
  int remainder = total_steps_0 % num_procs;

  int start_idx = rank * steps_per_proc + std::min(rank, remainder);
  int end_idx = start_idx + steps_per_proc + (rank < remainder ? 1 : 0) - 1;
  if (start_idx > total_steps_0) {
    start_idx = total_steps_0;
  }
  if (end_idx > total_steps_0) {
    end_idx = total_steps_0;
  }
  if (end_idx < start_idx) {
    end_idx = start_idx - 1;
  }

  double local_res = 0.0;

  int my_indices_count = end_idx - start_idx + 1;
  if (my_indices_count > 0) {
    int indices_per_thread = my_indices_count / num_omp_threads;
    int remainder_threads = my_indices_count % num_omp_threads;

#pragma omp parallel num_threads(num_omp_threads) default(none)                                                      \
    shared(start_idx, end_idx, indices_per_thread, remainder_threads, limits, actual_step, simpson_factor, function, \
               total_steps_0) reduction(+ : local_res)
    {
      int tid = omp_get_thread_num();
      int thread_start = start_idx + tid * indices_per_thread + std::min(tid, remainder_threads);
      int thread_end = thread_start + indices_per_thread + (tid < remainder_threads ? 1 : 0) - 1;
      if (thread_end > end_idx) {
        thread_end = end_idx;
      }

      for (int idx = thread_start; idx <= thread_end; ++idx) {
        double x0 = limits[0].first + idx * actual_step[0];
        double weight0 = Weight(idx, total_steps_0);

        local_res += IntegrateRemainingDims(x0, simpson_factor, limits, actual_step, function) * weight0;
      }
    }
  }

  double global_res = 0.0;
  MPI_Reduce(&local_res, &global_res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    I_res_ = global_res;
  }

  return true;
}

bool VinyaikinaEMultidimIntegrSimpsonALL::PostProcessingImpl() {
  GetOutput() = I_res_;
  return true;
}

}  // namespace vinyaikina_e_multidimensional_integrals_simpson_method
