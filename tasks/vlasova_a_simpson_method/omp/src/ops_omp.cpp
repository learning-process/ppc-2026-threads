#include "vlasova_a_simpson_method/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "vlasova_a_simpson_method/common/include/common.hpp"

namespace vlasova_a_simpson_method {

VlasovaASimpsonMethodOMP::VlasovaASimpsonMethodOMP(InType in) : task_data_(std::move(in)) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetOutput() = 0.0;
}

bool VlasovaASimpsonMethodOMP::ValidationImpl() {
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

bool VlasovaASimpsonMethodOMP::PreProcessingImpl() {
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

void VlasovaASimpsonMethodOMP::ComputeWeight(const std::vector<int> &index, double &weight) const {
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

void VlasovaASimpsonMethodOMP::ComputePoint(const std::vector<int> &index, std::vector<double> &point) const {
  size_t dim = index.size();
  point.resize(dim);

  for (size_t i = 0; i < dim; ++i) {
    point[i] = task_data_.a[i] + (index[i] * h_[i]);
  }
}

bool VlasovaASimpsonMethodOMP::RunImpl() {
  size_t dim = task_data_.a.size();

  size_t total_points = 1;
  for (size_t i = 0; i < dim; ++i) {
    total_points *= static_cast<size_t>(dimensions_[i]);
  }

  double sum = 0.0;

#pragma omp parallel default(none) shared(dim, total_points) reduction(+ : sum)
  {
    std::vector<int> cur_index(dim, 0);
    std::vector<double> cur_point;
    double local_weight = 0.0;

#pragma omp for schedule(static)
    for (size_t idx = 0; idx < total_points; ++idx) {
      size_t temp_idx = idx;
      for (size_t i = 0; i < dim; ++i) {
        cur_index[i] = static_cast<int>(temp_idx % static_cast<int>(dimensions_[i]));
        temp_idx /= static_cast<int>(dimensions_[i]);
      }

      ComputeWeight(cur_index, local_weight);
      ComputePoint(cur_index, cur_point);
      sum += local_weight * task_data_.func(cur_point);
    }
  }

  double factor = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    factor *= h_[i] / 3.0;
  }

  result_ = sum * factor;
  GetOutput() = result_;

  return true;
}

bool VlasovaASimpsonMethodOMP::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace vlasova_a_simpson_method
