#include "ovsyannikov_n_simpson_method_stl/stl/include/ops_stl.hpp"

#include <cmath>
#include <functional>
#include <numeric>
#include <thread>
#include <vector>

#include "ovsyannikov_n_simpson_method_stl/common/include/common.hpp"
#include "util/include/util.hpp"

namespace ovsyannikov_n_simpson_method_stl {

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
    partial_results[thread_idx] = local_sum;
  };

  int chunk_size = (nx_l + 1) / num_threads;
  for (unsigned int t = 0; t < num_threads; ++t) {
    int start = t * chunk_size;
    int end = (t == num_threads - 1) ? (nx_l + 1) : (start + chunk_size);
    threads.emplace_back(worker, start, end, t);
  }

  for (auto &th : threads) {
    if (th.joinable()) {
      th.join();
    }
  }

  double total_sum = std::accumulate(partial_results.begin(), partial_results.end(), 0.0);

  res_ = (hx * hy / 9.0) * total_sum;
  return true;
}

bool OvsyannikovNSimpsonMethodSTL::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

}  // namespace ovsyannikov_n_simpson_method_stl
