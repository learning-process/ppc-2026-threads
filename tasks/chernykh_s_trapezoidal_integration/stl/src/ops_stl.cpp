#include "chernykh_s_trapezoidal_integration/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <thread>
#include <utility>
#include <vector>

#include "chernykh_s_trapezoidal_integration/common/include/common.hpp"
#include "util/include/util.hpp"

namespace chernykh_s_trapezoidal_integration {

ChernykhSTrapezoidalIntegrationSTL::ChernykhSTrapezoidalIntegrationSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ChernykhSTrapezoidalIntegrationSTL::ValidationImpl() {
  const auto &input = this->GetInput();
  if (input.limits.empty() || input.limits.size() != input.steps.size()) {
    return false;
  }
  return std::ranges::all_of(input.steps, [](int s) { return s > 0; });
}

bool ChernykhSTrapezoidalIntegrationSTL::PreProcessingImpl() {
  return true;
}

void ChernykhSTrapezoidalIntegrationSTL::MemberWorker(int64_t start, int64_t end, double &local_result,
                                                      const std::vector<double> &h, int dims) {
  const auto &input = GetInput();
  double sum = 0.0;
  std::vector<double> point(dims);

  for (int64_t i = start; i < end; i++) {
    int64_t local_index = i;
    double weight = 1.0;

    for (int j = dims - 1; j >= 0; j--) {
      int64_t idx = local_index % (input.steps[j] + 1);
      local_index /= (input.steps[j] + 1);

      point[j] = input.limits[j].first + (static_cast<double>(idx) * h[j]);

      if (idx == 0 || idx == input.steps[j]) {
        weight *= 0.5;
      }
    }
    sum += input.func(point) * weight;
  }
  local_result = sum;
}

bool ChernykhSTrapezoidalIntegrationSTL::RunImpl() {
  auto &input = GetInput();
  int dims = static_cast<int>(input.limits.size());
  int64_t total_point = 1;
  for (int setka : input.steps) {
    total_point *= static_cast<int64_t>(setka) + 1;
  }

  std::vector<double> h(dims);
  for (int i = 0; i < dims; ++i) {
    h[i] = (input.limits[i].second - input.limits[i].first) / static_cast<double>(input.steps[i]);
  }

  int num_threads = ppc::util::GetNumThreads();
  if (num_threads == 0) {
    num_threads = 1;
  }

  std::vector<double> result(num_threads, 0.0);
  std::vector<std::thread> threads(num_threads);
  std::vector<std::pair<int64_t, int64_t>> borders(num_threads);

  int64_t points_per_thread = total_point / num_threads;
  int64_t remainder = total_point % num_threads;

  int64_t start = 0;
  for (int i = 0; i < num_threads; i++) {
    borders[i].first = start;
    start += points_per_thread;
    if (i < remainder) {
      start++;
    }
    borders[i].second = start;
  }

  auto work_function = [&](int64_t s, int64_t e, double &res) { MemberWorker(s, e, res, h, dims); };

  for (int i = 0; i < num_threads; ++i) {
    threads[i] = std::thread(work_function, borders[i].first, borders[i].second, std::ref(result[i]));
  }

  for (int i = 0; i < num_threads; ++i) {
    threads[i].join();
  }

  double output_result = 0.0;
  for (double loc_res : result) {
    output_result += loc_res;
  }

  for (int i = 0; i < dims; i++) {
    output_result *= h[i];
  }

  GetOutput() = output_result;
  return true;
}

bool ChernykhSTrapezoidalIntegrationSTL::PostProcessingImpl() {
  return true;
}

}  // namespace chernykh_s_trapezoidal_integration
