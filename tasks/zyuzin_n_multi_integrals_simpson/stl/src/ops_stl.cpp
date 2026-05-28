#include "zyuzin_n_multi_integrals_simpson/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <thread>
#include <vector>

#include "util/include/util.hpp"
#include "zyuzin_n_multi_integrals_simpson/common/include/common.hpp"

namespace zyuzin_n_multi_integrals_simpson {

ZyuzinNSimpsonSTL::ZyuzinNSimpsonSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ZyuzinNSimpsonSTL::ValidationImpl() {
  const auto &input = GetInput();
  if (input.lower_bounds.size() != input.upper_bounds.size() || input.lower_bounds.size() != input.n_steps.size()) {
    return false;
  }
  if (input.lower_bounds.empty()) {
    return false;
  }
  for (size_t i = 0; i < input.lower_bounds.size(); ++i) {
    if (input.lower_bounds[i] > input.upper_bounds[i]) {
      return false;
    }
    if (input.n_steps[i] <= 0 || input.n_steps[i] % 2 != 0) {
      return false;
    }
  }
  return static_cast<bool>(input.func);
}

bool ZyuzinNSimpsonSTL::PreProcessingImpl() {
  result_ = 0.0;
  return true;
}

double ZyuzinNSimpsonSTL::GetSimpsonWeight(int index, int n) {
  if (index == 0 || index == n) {
    return 1.0;
  }
  if (index % 2 == 1) {
    return 4.0;
  }
  return 2.0;
}

double ZyuzinNSimpsonSTL::ComputeSimpsonMultiDim() {
  const auto &input = GetInput();
  const size_t num_dims = input.lower_bounds.size();

  std::vector<double> h(num_dims);
  for (size_t dim = 0; dim < num_dims; ++dim) {
    h[dim] = (input.upper_bounds[dim] - input.lower_bounds[dim]) / input.n_steps[dim];
  }

  size_t total_points = 1;
  for (size_t dim = 0; dim < num_dims; ++dim) {
    total_points *= static_cast<size_t>(input.n_steps[dim]) + 1U;
  }

  const size_t requested_threads = static_cast<size_t>(std::max(1, ppc::util::GetNumThreads()));
  const size_t num_threads = std::min(requested_threads, total_points);

  std::vector<double> partial_sums(num_threads, 0.0);
  std::vector<std::thread> workers;
  workers.reserve(num_threads > 0 ? num_threads - 1 : 0);

  auto compute_chunk = [&](size_t thread_id) {
    const size_t begin = (thread_id * total_points) / num_threads;
    const size_t end = ((thread_id + 1) * total_points) / num_threads;

    std::vector<double> point(num_dims);
    double local_sum = 0.0;

    for (size_t point_idx = begin; point_idx < end; ++point_idx) {
      auto temp = point_idx;
      double weight = 1.0;

      for (size_t dim = 0; dim < num_dims; ++dim) {
        const auto axis_points = static_cast<size_t>(input.n_steps[dim]) + 1U;
        const auto index = static_cast<int>(temp % axis_points);
        temp /= axis_points;
        point[dim] = input.lower_bounds[dim] + (static_cast<double>(index) * h[dim]);
        weight *= GetSimpsonWeight(index, input.n_steps[dim]);
      }

      local_sum += weight * input.func(point);
    }

    partial_sums[thread_id] = local_sum;
  };

  for (size_t thread_id = 1; thread_id < num_threads; ++thread_id) {
    workers.emplace_back(compute_chunk, thread_id);
  }
  compute_chunk(0);

  for (auto &worker : workers) {
    worker.join();
  }

  const double sum = std::accumulate(partial_sums.begin(), partial_sums.end(), 0.0);

  double factor = 1.0;
  for (size_t dim = 0; dim < num_dims; ++dim) {
    factor *= h[dim] / 3.0;
  }

  return sum * factor;
}

bool ZyuzinNSimpsonSTL::RunImpl() {
  result_ = ComputeSimpsonMultiDim();
  return true;
}

bool ZyuzinNSimpsonSTL::PostProcessingImpl() {
  GetOutput() = result_;
  return true;
}

}  // namespace zyuzin_n_multi_integrals_simpson
