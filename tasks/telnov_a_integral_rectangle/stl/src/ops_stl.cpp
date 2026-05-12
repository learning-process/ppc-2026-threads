#include "telnov_a_integral_rectangle/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <thread>
#include <vector>

#include "telnov_a_integral_rectangle/common/include/common.hpp"
#include "util/include/util.hpp"

namespace telnov_a_integral_rectangle {

TelnovAIntegralRectangleSTL::TelnovAIntegralRectangleSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool TelnovAIntegralRectangleSTL::ValidationImpl() {
  return GetInput().first > 0 && GetInput().second > 0;
}

bool TelnovAIntegralRectangleSTL::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool TelnovAIntegralRectangleSTL::RunImpl() {
  const int n = GetInput().first;
  const int d = GetInput().second;

  const double a = 0.0;
  const double b = 1.0;
  const double h = (b - a) / static_cast<double>(n);

  const auto total_points = static_cast<int64_t>(std::pow(n, d));

  int thread_count = ppc::util::GetNumThreads();
  thread_count = std::max(1, std::min(thread_count, static_cast<int>(total_points)));

  std::vector<std::thread> threads(thread_count);
  std::vector<double> partial_sums(thread_count, 0.0);

  const int64_t block = total_points / thread_count;
  const int64_t remainder = total_points % thread_count;

  auto calculate_part = [n, d, a, h](int64_t begin, int64_t end) {
    double local_sum = 0.0;

    for (int64_t idx = begin; idx < end; ++idx) {
      int64_t current = idx;
      double f_value = 0.0;

      for (int dim = 0; dim < d; ++dim) {
        const int coord_index = static_cast<int>(current % n);
        current /= n;

        const double x = a + ((static_cast<double>(coord_index) + 0.5) * h);
        f_value += x;
      }

      local_sum += f_value;
    }

    return local_sum;
  };

  int64_t begin = 0;
  for (int i = 0; i < thread_count; ++i) {
    const int64_t current_block = block + (i < remainder ? 1 : 0);
    const int64_t end = begin + current_block;

    threads[i] = std::thread([&, i, begin, end]() { partial_sums[i] = calculate_part(begin, end); });

    begin = end;
  }

  for (auto &thread : threads) {
    thread.join();
  }

  double result = 0.0;
  for (const auto &value : partial_sums) {
    result += value;
  }

  GetOutput() = result * std::pow(h, d);
  return true;
}

bool TelnovAIntegralRectangleSTL::PostProcessingImpl() {
  return true;
}

}  // namespace telnov_a_integral_rectangle
