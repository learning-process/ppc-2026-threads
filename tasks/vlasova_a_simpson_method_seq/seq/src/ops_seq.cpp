#include "vlasova_a_simpson_method_seq/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "vlasova_a_simpson_method_seq/common/include/common.hpp"

namespace vlasova_a_simpson_method_seq {

VlasovaASimpsonMethodSEQ::VlasovaASimpsonMethodSEQ(const InType &in) : task_data_(in), result_(0.0) {
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

void VlasovaASimpsonMethodSEQ::NextIndex(std::vector<int> &Index) {
  size_t dim = Index.size();
  for (size_t i = 0; i < dim; ++i) {
    Index[i]++;
    if (Index[i] < dimensions_[i]) {
      return;
    }
    Index[i] = 0;
  }
}

double VlasovaASimpsonMethodSEQ::GetWeight(const std::vector<int> &Index) const {
  double weight = 1.0;
  size_t dim = Index.size();

  for (size_t i = 0; i < dim; ++i) {
    int idx = Index[i];
    int steps = task_data_.n[i];

    if (idx == 0 || idx == steps) {
      weight *= 1.0;
    } else if (idx % 2 == 0) {
      weight *= 2.0;
    } else {
      weight *= 4.0;
    }
  }

  return weight;
}

std::vector<double> VlasovaASimpsonMethodSEQ::GetPoint(const std::vector<int> &Index) const {
  size_t dim = Index.size();
  std::vector<double> point(dim);

  for (size_t i = 0; i < dim; ++i) {
    point[i] = task_data_.a[i] + (Index[i] * h_[i]);
  }

  return point;
}

bool VlasovaASimpsonMethodSEQ::RunImpl() {
  size_t dim = task_data_.a.size();
  std::vector<int> cur_Index(dim, 0);

  double sum = 0.0;
  bool has_more = true;

  while (has_more) {
    double weight = GetWeight(cur_Index);
    std::vector<double> cur_point = GetPoint(cur_Index);
    sum += weight * task_data_.func(cur_point);
    NextIndex(cur_Index);
    has_more = false;
    for (size_t i = 0; i < dim; ++i) {
      if (cur_Index[i] != 0) {
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
