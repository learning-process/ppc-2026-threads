#include "tsibareva_e_integral_calculate_trapezoid_method/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <thread>
#include <vector>

#include "tsibareva_e_integral_calculate_trapezoid_method/common/include/common.hpp"
#include "util/include/util.hpp"

namespace tsibareva_e_integral_calculate_trapezoid_method {

TsibarevaEIntegralCalculateTrapezoidMethodSTL::TsibarevaEIntegralCalculateTrapezoidMethodSTL(const InType &in)
    : BaseTask() {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodSTL::ValidationImpl() {
  return true;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodSTL::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodSTL::RunImpl() {
  int dim = GetInput().dim;

  std::vector<double> h(dim);
  std::vector<int> sizes(dim);
  int total_nodes = 1;
  for (int i = 0; i < dim; ++i) {
    h[i] = (GetInput().hi[i] - GetInput().lo[i]) / GetInput().steps[i];
    sizes[i] = GetInput().steps[i] + 1;
    total_nodes *= sizes[i];
  }

  const int num_threads = std::max<int>(1, std::min<int>(ppc::util::GetNumThreads(), total_nodes));
  std::vector<std::thread> threads(num_threads);
  std::vector<double> partial_sums(num_threads, 0.0);

  auto worker = [&](int thread_id, int start, int end) { MWork(thread_id, start, end, sizes, h, partial_sums, dim); };

  int nodes_per_thread = total_nodes / num_threads;
  int remainder_nodes = total_nodes % num_threads;
  int start = 0;
  for (int tid = 0; tid < num_threads; ++tid) {
    int end = start + nodes_per_thread + (tid < remainder_nodes ? 1 : 0);
    threads[tid] = std::thread(worker, tid, start, end);
    start = end;
  }

  for (auto &t : threads) {
    t.join();
  }

  double global_sum = 0.0;
  for (double s : partial_sums) {
    global_sum += s;
  }

  double res_h = 1.0;
  for (int i = 0; i < dim; ++i) {
    res_h *= h[i];
  }

  GetOutput() = global_sum * res_h;
  return true;
}

bool TsibarevaEIntegralCalculateTrapezoidMethodSTL::PostProcessingImpl() {
  return true;
}

void TsibarevaEIntegralCalculateTrapezoidMethodSTL::MWork(int thread_id, int start, int end,
                                                          const std::vector<int> &sizes, const std::vector<double> &h,
                                                          std::vector<double> &partial_sums, int dim) {
  double local_sum = 0.0;

  for (int node = start; node < end; ++node) {
    int remainder = node;
    double node_weight = 1.0;
    std::vector<double> point(dim);

    for (int i = dim - 1; i >= 0; --i) {
      int idx = remainder % sizes[i];
      remainder /= sizes[i];

      if (idx == 0 || idx == GetInput().steps[i]) {
        node_weight *= 0.5;
      }

      point[i] = GetInput().lo[i] + (idx * h[i]);
    }

    local_sum += node_weight * GetInput().f(point);
  }
  partial_sums[thread_id] = local_sum;
}

}  // namespace tsibareva_e_integral_calculate_trapezoid_method
