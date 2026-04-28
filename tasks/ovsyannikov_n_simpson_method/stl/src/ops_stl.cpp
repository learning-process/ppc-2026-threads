#include "ovsyannikov_n_simpson_method/stl/include/ops_stl.hpp"

#include <cmath>
#include <cstddef>
#include <numeric>
#include <thread>
#include <utility>
#include <vector>

#include "ovsyannikov_n_simpson_method/common/include/common.hpp"
#include "util/include/util.hpp"

namespace ovsyannikov_n_simpson_method {

double OvsyannikovNSimpsonMethodSTL::Function(double x, double y) {
  return x + y;
}

double OvsyannikovNSimpsonMethodSTL::GetCoeff(int i, int n) {
  if (i == 0 || i == n) {
    return 1.0;
  }
  return (i % 2 == 1) ? 4.0 : 2.0;
}

OvsyannikovNSimpsonMethodSTL::OvsyannikovNSimpsonMethodSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool OvsyannikovNSimpsonMethodSTL::ValidationImpl() {
  return GetInput().nx > 0 && GetInput().nx % 2 == 0 && GetInput().ny > 0 && GetInput().ny % 2 == 0;
}

bool OvsyannikovNSimpsonMethodSTL::PreProcessingImpl() {
  params_ = GetInput();
  res_ = 0.0;
  return true;
}

bool OvsyannikovNSimpsonMethodSTL::RunImpl() {
  const int nx_l = params_.nx;
  const int ny_l = params_.ny;
  const double ax_l = params_.ax;
  const double ay_l = params_.ay;
  const double hx = (params_.bx - params_.ax) / nx_l;
  const double hy = (params_.by - params_.ay) / ny_l;

  unsigned int num_threads = ppc::util::GetNumThreads();
  if (num_threads == 0) {
    num_threads = 2;
  }

  std::vector<double> partial_results(num_threads, 0.0);
  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  auto worker = [&](int start_i, int end_i, int thread_idx) {
    double local_sum = 0.0;
    for (int i = start_i; i < end_i; ++i) {
      const double x = ax_l + (static_cast<double>(i) * hx);
      const double coeff_x = GetCoeff(i, nx_l);
      double row_sum = 0.0;
      for (int j = 0; j <= ny_l; ++j) {
        const double y = ay_l + (static_cast<double>(j) * hy);
        const double coeff_y = GetCoeff(j, ny_l);
        row_sum += coeff_y * Function(x, y);
      }
      local_sum += coeff_x * row_sum;
    }
    partial_results[static_cast<std::size_t>(thread_idx)] = local_sum;
  };

  const int total_tasks = nx_l + 1;
  const int chunk_size = total_tasks / static_cast<int>(num_threads);
  const int extra = total_tasks % static_cast<int>(num_threads);

  int current_start = 0;
  for (unsigned int i_thread = 0; i_thread < num_threads; ++i_thread) {
    int current_end = current_start + chunk_size + (std::cmp_less(i_thread, extra) ? 1 : 0);
    threads.emplace_back(worker, current_start, current_end, static_cast<int>(i_thread));
    current_start = current_end;
  }

  for (auto &th : threads) {
    if (th.joinable()) {
      th.join();
    }
  }

  res_ = (hx * hy / 9.0) * std::accumulate(partial_results.begin(), partial_results.end(), 0.0);
  return true;
}

bool OvsyannikovNSimpsonMethodSTL::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

}  // namespace ovsyannikov_n_simpson_method
