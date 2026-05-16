#include "vlasova_a_simpson_method/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <execution>
#include <functional>
#include <numeric>
#include <ranges>
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

double VlasovaASimpsonMethodSTL::TraverseFrom(size_t dim, std::vector<int> &cur_index, std::vector<double> &cur_point,
                                              size_t d) const {
  if (d == dim) {
    double weight = 0.0;
    ComputeWeight(cur_index, weight);
    ComputePoint(cur_index, cur_point);
    return weight * task_data_.func(cur_point);
  }
  double local_sum = 0.0;
  for (int i = 0; i < dimensions_[d]; ++i) {
    cur_index[d] = i;
    local_sum += TraverseFrom(dim, cur_index, cur_point, d + 1);
  }
  return local_sum;
}

bool VlasovaASimpsonMethodSTL::RunImpl() {
  size_t dim = task_data_.a.size();
  int first_dim_size = dimensions_[0];

  std::vector<int> indices(first_dim_size);
  std::ranges::iota(indices, 0);

  std::vector<double> partial_sums(first_dim_size, 0.0);

  std::transform(std::execution::par, indices.begin(), indices.end(), partial_sums.begin(), [this, dim](int idx0) {
    thread_local std::vector<int> cur_index;
    thread_local std::vector<double> cur_point;
    cur_index.assign(dim, 0);
    cur_index[0] = idx0;
    return TraverseFrom(dim, cur_index, cur_point, 1);
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
