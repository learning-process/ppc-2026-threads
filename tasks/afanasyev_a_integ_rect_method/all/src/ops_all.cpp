#include "afanasyev_a_integ_rect_method/all/include/ops_all.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <numeric>
#include <thread>
#include <vector>

#include "afanasyev_a_integ_rect_method/common/include/common.hpp"
#include "util/include/util.hpp"

namespace afanasyev_a_integ_rect_method {
namespace {

double ExampleIntegrand(const std::array<double, 3> &x) {
  double s = 0.0;
  for (double xi : x) {
    s += xi * xi;
  }
  return std::exp(-s);
}

double ComputeSumForI(int i, int n, double h) {
  double sum = 0.0;
  std::array<double, 3> x{};
  x[0] = (static_cast<double>(i) + 0.5) * h;
  for (int j = 0; j < n; ++j) {
    x[1] = (static_cast<double>(j) + 0.5) * h;
    for (int k = 0; k < n; ++k) {
      x[2] = (static_cast<double>(k) + 0.5) * h;
      sum += ExampleIntegrand(x);
    }
  }
  return sum;
}

}  // namespace

AfanasyevAIntegRectMethodALL::AfanasyevAIntegRectMethodALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool AfanasyevAIntegRectMethodALL::ValidationImpl() {
  return (GetInput() > 0);
}

bool AfanasyevAIntegRectMethodALL::PreProcessingImpl() {
  return true;
}

bool AfanasyevAIntegRectMethodALL::RunImpl() {
  const int n = GetInput();
  if (n <= 0) {
    return false;
  }

  const int k_dim = 3;
  const double h = 1.0 / static_cast<double>(n);
  const double volume = std::pow(h, k_dim);

  const int num_threads = std::max(1, ppc::util::GetNumThreads());
  const int worker_count = std::min(num_threads, n);

  std::vector<std::thread> workers;
  std::vector<double> local_sums(worker_count, 0.0);

  const int chunk_size = n / worker_count;
  const int remainder = n % worker_count;
  int start_i = 0;

  for (int thread_id = 0; thread_id < worker_count; ++thread_id) {
    const int end_i = start_i + chunk_size + (thread_id < remainder ? 1 : 0);
    workers.emplace_back([thread_id, start_i, end_i, n, h, &local_sums]() {
      local_sums[thread_id] = tbb::parallel_reduce(tbb::blocked_range<int>(start_i, end_i), 0.0,
                                                   [&](const tbb::blocked_range<int> &range, double sum) {
        for (int i = range.begin(); i < range.end(); ++i) {
          sum += ComputeSumForI(i, n, h);
        }
        return sum;
      }, std::plus<>());
    });
    start_i = end_i;
  }

  for (auto &worker : workers) {
    worker.join();
  }

  const double sum = std::accumulate(local_sums.begin(), local_sums.end(), 0.0);
  GetOutput() = sum * volume;
  return true;
}

bool AfanasyevAIntegRectMethodALL::PostProcessingImpl() {
  return true;
}

}  // namespace afanasyev_a_integ_rect_method
