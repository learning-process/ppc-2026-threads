#include "vlasova_a_simpson_method/stl/include/ops_stl.hpp"

#include <cmath>
#include <cstddef>
#include <numeric>
#include <thread>
#include <utility>
#include <vector>

#include "util/include/util.hpp"
#include "vlasova_a_simpson_method/common/include/common.hpp"

namespace vlasova_a_simpson_method {

VlasovaASimpsonMethodSTL::VlasovaASimpsonMethodSTL(InType in) : task_data_(std::move(in)) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetOutput() = 0.0;
}

bool VlasovaASimpsonMethodSTL::ValidationImpl() {
  size_t dim = task_data_.a.size();

  if (dim == 0 || dim != task_data_.b.size() || dim != task_data_.n.size()) {
    return false;
  }

  for (size_t i = 0; i < dim; ++i) {
    if (task_data_.a[i] >= task_data_.b[i]) {
      return false;
    }
    if (task_data_.n[i] <= 0 || task_data_.n[i] % 2 != 0) {
      return false;
    }
  }

  if (!task_data_.func) {
    return false;
  }

  return GetOutput() == 0.0;
}

bool VlasovaASimpsonMethodSTL::PreProcessingImpl() {
  result_ = 0.0;
  GetOutput() = 0.0;

  size_t dim = task_data_.a.size();
  h_.resize(dim);
  dimensions_.resize(dim);

  for (size_t i = 0; i < dim; ++i) {
    h_[i] = (task_data_.b[i] - task_data_.a[i]) / task_data_.n[i];
    dimensions_[i] = task_data_.n[i] + 1;
  }

  return true;
}

void VlasovaASimpsonMethodSTL::ComputeWeight(const std::vector<int> &index, double &weight) const {
  weight = 1.0;
  size_t dim = index.size();

  for (size_t i = 0; i < dim; ++i) {
    int idx = index[i];
    int steps = task_data_.n[i];

    if (idx == 0 || idx == steps) {
      weight *= 1.0;
    } else if (idx % 2 == 0) {
      weight *= 2.0;
    } else {
      weight *= 4.0;
    }
  }
}

void VlasovaASimpsonMethodSTL::ComputePoint(const std::vector<int> &index, std::vector<double> &point) const {
  size_t dim = index.size();
  point.resize(dim);

  for (size_t i = 0; i < dim; ++i) {
    point[i] = task_data_.a[i] + (index[i] * h_[i]);
  }
}

void VlasovaASimpsonMethodSTL::ComputePartialSumRange(int start_idx, int end_idx, double &partial_sum) const {
  size_t dim = task_data_.a.size();
  std::vector<int> cur_index(dim, 0);
  std::vector<double> cur_point;
  double local_sum = 0.0;

  for (int idx = start_idx; idx < end_idx; ++idx) {
    auto temp_idx = static_cast<size_t>(idx);
    for (size_t i = 0; i < dim; ++i) {
      cur_index[i] = static_cast<int>(temp_idx % static_cast<size_t>(dimensions_[i]));
      temp_idx /= static_cast<size_t>(dimensions_[i]);
    }

    double weight = 0.0;
    ComputeWeight(cur_index, weight);
    ComputePoint(cur_index, cur_point);
    local_sum += weight * task_data_.func(cur_point);
  }

  partial_sum = local_sum;
}

void VlasovaASimpsonMethodSTL::ProcessSequential(double &sum) const {
  size_t dim = task_data_.a.size();
  std::vector<int> cur_index(dim, 0);
  std::vector<double> cur_point;
  sum = 0.0;

  size_t total_points = 1;
  for (size_t i = 0; i < dim; ++i) {
    total_points *= static_cast<size_t>(dimensions_[i]);
  }

  for (size_t idx = 0; idx < total_points; ++idx) {
    auto temp_idx = idx;
    for (size_t i = 0; i < dim; ++i) {
      cur_index[i] = static_cast<int>(temp_idx % static_cast<size_t>(dimensions_[i]));
      temp_idx /= static_cast<size_t>(dimensions_[i]);
    }

    double weight = 0.0;
    ComputeWeight(cur_index, weight);
    ComputePoint(cur_index, cur_point);
    sum += weight * task_data_.func(cur_point);
  }
}

void VlasovaASimpsonMethodSTL::ProcessParallel(double &sum) const {
  size_t dim = task_data_.a.size();

  size_t total_points = 1;
  for (size_t i = 0; i < dim; ++i) {
    total_points *= static_cast<size_t>(dimensions_[i]);
  }

  unsigned int num_threads = ppc::util::GetNumThreads();
  if (num_threads == 0) {
    num_threads = std::thread::hardware_concurrency();
  }
  if (num_threads == 0) {
    num_threads = 2;
  }

  if (static_cast<size_t>(num_threads) > total_points) {
    num_threads = static_cast<unsigned int>(total_points);
  }

  size_t points_per_thread = total_points / num_threads;
  size_t remainder = total_points % num_threads;

  std::vector<std::thread> threads;
  std::vector<double> partial_sums(num_threads, 0.0);

  size_t start_idx = 0;
  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    size_t end_idx = start_idx + points_per_thread + (thread_idx < remainder ? 1 : 0);
    threads.emplace_back([this, start_idx, end_idx, &partial_sums, thread_idx]() {
      double local_sum = 0.0;
      ComputePartialSumRange(static_cast<int>(start_idx), static_cast<int>(end_idx), local_sum);
      partial_sums[thread_idx] = local_sum;
    });
    start_idx = end_idx;
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }

  sum = std::accumulate(partial_sums.begin(), partial_sums.end(), 0.0);
}

bool VlasovaASimpsonMethodSTL::RunImpl() {
  size_t dim = task_data_.a.size();

  size_t total_points = 1;
  for (size_t i = 0; i < dim; ++i) {
    total_points *= static_cast<size_t>(dimensions_[i]);
  }

  double sum = 0.0;
  if (total_points < 10000) {
    ProcessSequential(sum);
  } else {
    ProcessParallel(sum);
  }

  double factor = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    factor *= h_[i] / 3.0;
  }

  result_ = sum * factor;
  GetOutput() = result_;

  return true;
}

bool VlasovaASimpsonMethodSTL::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace vlasova_a_simpson_method
