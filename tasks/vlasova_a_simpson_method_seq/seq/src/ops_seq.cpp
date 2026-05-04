#include "vlasova_a_simpson_method_seq/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

#include "vlasova_a_simpson_method_seq/common/include/common.hpp"

namespace vlasova_a_simpson_method_seq {

VlasovaASimpsonMethodSEQ::VlasovaASimpsonMethodSEQ(InType in) : task_data_(std::move(in)) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetOutput() = 0.0;
}

bool VlasovaASimpsonMethodSEQ::ValidationImpl() {
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

bool VlasovaASimpsonMethodSEQ::PreProcessingImpl() {
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

void VlasovaASimpsonMethodSEQ::Nextindex(std::vector<int> &index) {
  size_t dim = index.size();
  for (size_t i = 0; i < dim; ++i) {
    index[i]++;
    if (index[i] < dimensions_[i]) {
      return;
    }
    index[i] = 0;
  }
}

void VlasovaASimpsonMethodSEQ::ComputeWeight(const std::vector<int> &index, double &weight) const {
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

void VlasovaASimpsonMethodSEQ::ComputePoint(const std::vector<int> &index, std::vector<double> &point) const {
  size_t dim = index.size();
  point.resize(dim);

  for (size_t i = 0; i < dim; ++i) {
    point[i] = task_data_.a[i] + (index[i] * h_[i]);
  }
}

bool VlasovaASimpsonMethodSEQ::RunImpl() {
  size_t dim = task_data_.a.size();

  std::vector<int> cur_index(dim, 0);
  std::vector<double> cur_point;
  double sum = 0.0;
  bool has_more = true;

  while (has_more) {
    double weight = 0.0;

    ComputeWeight(cur_index, weight);
    ComputePoint(cur_index, cur_point);

    sum += weight * task_data_.func(cur_point);
    Nextindex(cur_index);
    has_more = false;
    for (size_t i = 0; i < dim; ++i) {
      if (cur_index[i] != 0) {
        has_more = true;
        break;
      }
    }
  }

  // Множитель: (h1 * h2 * ... * hd) / 3^d
  double factor = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    factor *= h_[i] / 3.0;
  }

  result_ = sum * factor;
  GetOutput() = result_;

  return true;
}

bool VlasovaASimpsonMethodSEQ::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace vlasova_a_simpson_method_seq
