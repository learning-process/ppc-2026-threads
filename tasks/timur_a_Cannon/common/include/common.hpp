#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace timur_a_Cannon {

struct Matrix {
  std::size_t n{};
  std::vector<double> data;

  Matrix() = default;
  Matrix(std::size_t size, double init_val = 0.0) : n(size), data(size * size, init_val) {}

  double &operator()(std::size_t row, std::size_t col) {
    return data[row * n + col];
  }
  const double &operator()(std::size_t row, std::size_t col) const {
    return data[row * n + col];
  }

  bool operator==(const Matrix &other) const {
    return n == other.n && data == other.data;
  }
  bool operator!=(const Matrix &other) const {
    return !(*this == other);
  }
};

struct TaskData {
  Matrix A;
  Matrix B;

  bool operator==(const TaskData &other) const {
    return A == other.A && B == other.B;
  }
  bool operator!=(const TaskData &other) const {
    return !(*this == other);
  }
};

using InType = TaskData;
using OutType = Matrix;
using TestType = std::tuple<std::size_t, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace timur_a_Cannon
