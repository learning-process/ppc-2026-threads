#include "vlasova_a_simpson_method_seq/seq/include/ops_seq.hpp"

#include <cmath>
#include <functional>
#include <stdexcept>

#include "util/include/util.hpp"
#include "vlasova_a_simpson_method_seq/common/include/common.hpp"

namespace vlasova_a_simpson_method_seq {

VlasovaASimpsonMethodSEQ::VlasovaASimpsonMethodSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  task_data_ = in;
  result_ = 0.0;
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
  for (size_t i = 0; i < dim; ++i) {
    h_[i] = (task_data_.b[i] - task_data_.a[i]) / task_data_.n[i];
  }

  return true;
}

double VlasovaASimpsonMethodSEQ::SimpsonRecursive(size_t dim, std::vector<double> &point) {
  size_t total_dims = task_data_.a.size();

  if (dim == total_dims) {
    return task_data_.func(point);
  }

  double result = 0.0;
  int steps = task_data_.n[dim];

  for (int i = 0; i <= steps; ++i) {
    point[dim] = task_data_.a[dim] + i * h_[dim];

    double weight = 1.0;
    if (i == 0 || i == steps) {
      weight = 1.0;
    } else if (i % 2 == 0) {
      weight = 2.0;
    } else {
      weight = 4.0;
    }

    result += weight * SimpsonRecursive(dim + 1, point);
  }

  return result;
}

bool VlasovaASimpsonMethodSEQ::RunImpl() {
  size_t dim = task_data_.a.size();
  std::vector<double> point(dim, 0.0);

  double sum = SimpsonRecursive(0, point);

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
