#include "zyuzin_n_multi_integrals_simpson/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <cstddef>
#include <thread>
#include <vector>

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/global_control.h"
#include "oneapi/tbb/parallel_reduce.h"
#include "zyuzin_n_multi_integrals_simpson/common/include/common.hpp"

namespace zyuzin_n_multi_integrals_simpson {

ZyuzinNSimpsonTBB::ZyuzinNSimpsonTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ZyuzinNSimpsonTBB::ValidationImpl() {
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

bool ZyuzinNSimpsonTBB::PreProcessingImpl() {
  result_ = 0.0;
  return true;
}

double ZyuzinNSimpsonTBB::GetSimpsonWeight(int index, int n) {
  if (index == 0 || index == n) {
    return 1.0;
  }
  if (index % 2 == 1) {
    return 4.0;
  }
  return 2.0;
}

double ZyuzinNSimpsonTBB::ComputeSimpsonMultiDim() {
  const auto &input = GetInput();
  const size_t num_dims = input.lower_bounds.size();

  std::vector<double> h(num_dims);
  for (size_t dim = 0; dim < num_dims; ++dim) {
    h[dim] = (input.upper_bounds[dim] - input.lower_bounds[dim]) / input.n_steps[dim];
  }

  size_t total_points = 1;
  for (size_t dim = 0; dim < num_dims; ++dim) {
    total_points *= static_cast<size_t>(input.n_steps[dim] + 1);
  }

  const int max_threads = static_cast<int>(std::max(1U, std::thread::hardware_concurrency()));
  oneapi::tbb::global_control thread_limiter(oneapi::tbb::global_control::max_allowed_parallelism, max_threads);
  constexpr size_t kGrainSize = 4096;
  const double sum =
      oneapi::tbb::parallel_reduce(oneapi::tbb::blocked_range<size_t>(0, total_points, kGrainSize), 0.0,
                                   [&](const oneapi::tbb::blocked_range<size_t> &range, double local_sum) {
    std::vector<double> point(num_dims);
    for (size_t point_idx = range.begin(); point_idx < range.end(); ++point_idx) {
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
    return local_sum;
  }, [](double left, double right) { return left + right; });

  double factor = 1.0;
  for (size_t dim = 0; dim < num_dims; ++dim) {
    factor *= h[dim] / 3.0;
  }

  return sum * factor;
}

bool ZyuzinNSimpsonTBB::RunImpl() {
  result_ = ComputeSimpsonMultiDim();
  return true;
}

bool ZyuzinNSimpsonTBB::PostProcessingImpl() {
  GetOutput() = result_;
  return true;
}

}  // namespace zyuzin_n_multi_integrals_simpson
