#include "afanasyev_a_integ_rect_method/stl/include/ops_stl.hpp"

#include <cmath>
#include <functional>
#include <numeric>
#include <thread>
#include <vector>

#include "afanasyev_a_integ_rect_method/common/include/common.hpp"
#include "util/include/util.hpp"

namespace afanasyev_a_integ_rect_method {
namespace {

double ExampleIntegrand(const std::vector<double> &x) {
  double s = 0.0;
  for (double xi : x) {
    s += xi * xi;
  }
  return std::exp(-s);
}

void ComputePartialSum(int start, int end, int n, int k_dim, double h, double &local_sum) {
  std::vector<int> idx(k_dim, 0);
  std::vector<double> x(k_dim, 0.0);
  local_sum = 0.0;

  for (int i0 = start; i0 < end; ++i0) {
    idx[0] = i0;
    for (int i1 = 0; i1 < n; ++i1) {
      idx[1] = i1;
      for (int i2 = 0; i2 < n; ++i2) {
        idx[2] = i2;
        for (int dim = 0; dim < k_dim; ++dim) {
          x[dim] = (static_cast<double>(idx[dim]) + 0.5) * h;
        }
        local_sum += ExampleIntegrand(x);
      }
    }
  }
}

}  // namespace

AfanasyevAIntegRectMethodSTL::AfanasyevAIntegRectMethodSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool AfanasyevAIntegRectMethodSTL::ValidationImpl() {
  return (GetInput() > 0);
}

bool AfanasyevAIntegRectMethodSTL::PreProcessingImpl() {
  return true;
}

bool AfanasyevAIntegRectMethodSTL::RunImpl() {
  const int n = GetInput();
  if (n <= 0) {
    return false;
  }

  const int k_dim = 3;
  const double h = 1.0 / static_cast<double>(n);

  const int num_threads = ppc::util::GetNumThreads();
  std::vector<std::thread> threads;
  std::vector<double> partial_sums(num_threads, 0.0);

  int chunk_size = n / num_threads;
  int remainder = n % num_threads;

  int start = 0;
  for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
    int end = start + chunk_size + (thread_id < remainder ? 1 : 0);
    threads.emplace_back(ComputePartialSum, start, end, n, k_dim, h, std::ref(partial_sums[thread_id]));
    start = end;
  }

  for (auto &th : threads) {
    th.join();
  }

  double sum = std::accumulate(partial_sums.begin(), partial_sums.end(), 0.0);
  const double volume = std::pow(h, k_dim);
  GetOutput() = sum * volume;

  return true;
}

bool AfanasyevAIntegRectMethodSTL::PostProcessingImpl() {
  return true;
}

}  // namespace afanasyev_a_integ_rect_method
