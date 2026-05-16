#include "vlasova_a_simpson_method/stl/include/ops_stl.hpp"

#include <cmath>
#include <cstddef>
#include <numeric>
#include <thread>
#include <utility>
#include <vector>

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
    size_t temp_idx = static_cast<size_t>(idx);
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

bool VlasovaASimpsonMethodSTL::RunImpl() {
  size_t dim = task_data_.a.size();

  size_t total_points = 1;
  for (size_t i = 0; i < dim; ++i) {
    total_points *= static_cast<size_t>(dimensions_[i]);
  }

  // Для маленьких задач используем последовательное выполнение
  if (total_points < 10000) {
    std::vector<int> cur_index(dim, 0);
    std::vector<double> cur_point;
    double sum = 0.0;

    for (size_t idx = 0; idx < total_points; ++idx) {
      size_t temp_idx = idx;
      for (size_t i = 0; i < dim; ++i) {
        cur_index[i] = static_cast<int>(temp_idx % static_cast<size_t>(dimensions_[i]));
        temp_idx /= static_cast<size_t>(dimensions_[i]);
      }

      double weight = 0.0;
      ComputeWeight(cur_index, weight);
      ComputePoint(cur_index, cur_point);
      sum += weight * task_data_.func(cur_point);
    }

    double factor = 1.0;
    for (size_t i = 0; i < dim; ++i) {
      factor *= h_[i] / 3.0;
    }

    result_ = sum * factor;
    GetOutput() = result_;
    return true;
  }

  // Параллельное выполнение с std::thread (как у коллег)
  unsigned int num_threads = ppc::util::GetNumThreads();
  if (num_threads == 0) {
    num_threads = std::thread::hardware_concurrency();
  }
  if (num_threads == 0) {
    num_threads = 2;
  }

  // Не создаем больше потоков, чем нужно
  if (static_cast<size_t>(num_threads) > total_points) {
    num_threads = static_cast<unsigned int>(total_points);
  }

  size_t points_per_thread = total_points / num_threads;
  size_t remainder = total_points % num_threads;

  std::vector<std::thread> threads;
  std::vector<double> partial_sums(num_threads, 0.0);

  size_t start_idx = 0;
  for (unsigned int t = 0; t < num_threads; ++t) {
    size_t end_idx = start_idx + points_per_thread + (t < remainder ? 1 : 0);
    threads.emplace_back([this, start_idx, end_idx, &partial_sums, t]() {
      double local_sum = 0.0;
      ComputePartialSumRange(static_cast<int>(start_idx), static_cast<int>(end_idx), local_sum);
      partial_sums[t] = local_sum;
    });
    start_idx = end_idx;
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }

  double sum = std::accumulate(partial_sums.begin(), partial_sums.end(), 0.0);

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
