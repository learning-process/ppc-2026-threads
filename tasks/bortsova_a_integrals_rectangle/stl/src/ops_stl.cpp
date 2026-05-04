#include "bortsova_a_integrals_rectangle/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstdint>
#include <thread>
#include <vector>

#include "bortsova_a_integrals_rectangle/common/include/common.hpp"
#include "util/include/util.hpp"

namespace bortsova_a_integrals_rectangle {

BortsovaAIntegralsRectangleSTL::BortsovaAIntegralsRectangleSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool BortsovaAIntegralsRectangleSTL::ValidationImpl() {
  const auto &input = GetInput();
  return input.func && !input.lower_bounds.empty() && input.lower_bounds.size() == input.upper_bounds.size() &&
         input.num_steps > 0;
}

bool BortsovaAIntegralsRectangleSTL::PreProcessingImpl() {
  const auto &input = GetInput();
  func_ = input.func;
  num_steps_ = input.num_steps;
  dims_ = static_cast<int>(input.lower_bounds.size());

  midpoints_.resize(dims_);
  volume_ = 1.0;
  total_points_ = 1;

  for (int di = 0; di < dims_; di++) {
    double step = (input.upper_bounds[di] - input.lower_bounds[di]) / static_cast<double>(num_steps_);
    volume_ *= step;
    total_points_ *= num_steps_;

    midpoints_[di].resize(num_steps_);
    for (int si = 0; si < num_steps_; si++) {
      midpoints_[di][si] = input.lower_bounds[di] + ((si + 0.5) * step);
    }
  }

  return true;
}

double BortsovaAIntegralsRectangleSTL::ComputePartialSum(int64_t begin, int64_t end) {
  std::vector<int> indices(dims_, 0);
  std::vector<double> point(dims_);

  int64_t temp = begin;
  for (int di = dims_ - 1; di >= 0; di--) {
    indices[di] = static_cast<int>(temp % num_steps_);
    temp /= num_steps_;
  }

  double local_sum = 0.0;
  for (int64_t pt = begin; pt < end; pt++) {
    for (int di = 0; di < dims_; di++) {
      point[di] = midpoints_[di][indices[di]];
    }
    local_sum += func_(point);

    for (int di = dims_ - 1; di >= 0; di--) {
      indices[di]++;
      if (indices[di] < num_steps_) {
        break;
      }
      indices[di] = 0;
    }
  }
  return local_sum;
}

bool BortsovaAIntegralsRectangleSTL::RunImpl() {
  int num_threads = ppc::util::GetNumThreads();
  std::vector<double> partial_sums(num_threads, 0.0);

  int64_t chunk = total_points_ / num_threads;
  int64_t remainder = total_points_ % num_threads;

  std::vector<std::thread> threads(num_threads - 1);
  for (int ti = 1; ti < num_threads; ti++) {
    int64_t begin = (ti * chunk) + std::min(static_cast<int64_t>(ti), remainder);
    int64_t end = begin + chunk + (static_cast<int64_t>(ti) < remainder ? 1 : 0);
    threads[ti - 1] = std::thread([&, ti, begin, end]() { partial_sums[ti] = ComputePartialSum(begin, end); });
  }

  partial_sums[0] = ComputePartialSum(0, chunk + (remainder > 0 ? 1 : 0));

  double sum = partial_sums[0];
  for (int ti = 1; ti < num_threads; ti++) {
    threads[ti - 1].join();
    sum += partial_sums[ti];
  }

  GetOutput() = sum * volume_;
  return true;
}

bool BortsovaAIntegralsRectangleSTL::PostProcessingImpl() {
  return true;
}

}  // namespace bortsova_a_integrals_rectangle
