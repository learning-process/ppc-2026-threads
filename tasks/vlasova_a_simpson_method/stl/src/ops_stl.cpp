#include "vlasova_a_simpson_method/stl/include/ops_stl.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <execution>
#include <functional>
#include <numeric>
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

bool VlasovaASimpsonMethodSTL::RunImpl() {
  size_t dim = task_data_.a.size();

  size_t total_points = 1;
  for (size_t i = 0; i < dim; ++i) {
    total_points *= static_cast<size_t>(dimensions_[i]);
  }

  if (total_points < 10000) {
    std::vector<int> cur_index(dim, 0);
    std::vector<double> cur_point;
    double sum = 0.0;

    std::function<bool(size_t)> next_index = [&](size_t d) -> bool {
      if (d == dim) {
        return false;
      }
      cur_index[d]++;
      if (cur_index[d] < dimensions_[d]) {
        return true;
      }
      cur_index[d] = 0;
      return next_index(d + 1);
    };

    do {
      double weight = 0.0;
      ComputeWeight(cur_index, weight);
      ComputePoint(cur_index, cur_point);
      sum += weight * task_data_.func(cur_point);
    } while (next_index(0));

    double factor = 1.0;
    for (size_t i = 0; i < dim; ++i) {
      factor *= h_[i] / 3.0;
    }

    result_ = sum * factor;
    GetOutput() = result_;
    return true;
  }

  int first_dim_size = dimensions_[0];
  std::vector<double> partial_sums(first_dim_size);
  std::vector<int> indices(first_dim_size);
  std::iota(indices.begin(), indices.end(), 0);

  std::for_each(std::execution::par, indices.begin(), indices.end(), [this, dim, &partial_sums](int idx0) {
    static thread_local std::vector<int> cur_index(dim);
    static thread_local std::vector<double> cur_point;
    cur_index[0] = idx0;

    double local_sum = 0.0;

    // Рекурсивный обход остальных измерений
    std::function<void(size_t)> traverse = [&](size_t d) {
      if (d == dim) {
        double weight = 0.0;
        ComputeWeight(cur_index, weight);
        ComputePoint(cur_index, cur_point);
        local_sum += weight * task_data_.func(cur_point);
        return;
      }

      for (int i = 0; i < dimensions_[d]; ++i) {
        cur_index[d] = i;
        traverse(d + 1);
      }
    };

    traverse(1);  // Начинаем со второго измерения
    partial_sums[idx0] = local_sum;
  });

  double sum = std::reduce(std::execution::par_unseq, partial_sums.begin(), partial_sums.end(), 0.0);

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
