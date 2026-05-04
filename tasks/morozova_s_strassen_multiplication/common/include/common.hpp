#pragma once

#include <cmath>
#include <cstddef>  // для size_t
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace morozova_s_strassen_multiplication {
using InType = std::vector<double>;
using OutType = std::vector<double>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

struct Matrix {
  std::vector<double> data;
  int size;

  Matrix() : size(0) {}

  explicit Matrix(int n) : data(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0), size(n) {}

  double &operator()(int i, int j) {
    return data[(i * size) + j];
  }
  const double &operator()(int i, int j) const {
    return data[(i * size) + j];
  }

  bool operator==(const Matrix &other) const {
    if (size != other.size) {
      return false;
    }
    const double eps = 1e-6;
    for (size_t i = 0; i < data.size(); ++i) {
      if (std::abs(data[i] - other.data[i]) > eps) {
        return false;
      }
    }
    return true;
  }
};

}  // namespace morozova_s_strassen_multiplication
