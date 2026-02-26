#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace yakimov_i_mult_of_dense_matrices_fox_algorithm_seq {

using InType = int;
using OutType = double;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

struct DenseMatrix {
  std::vector<double> data;
  int rows = 0;
  int cols = 0;

  DenseMatrix() = default;

  double &operator()(int i, int j) {
    std::size_t index = (static_cast<std::size_t>(i) * static_cast<std::size_t>(cols)) + static_cast<std::size_t>(j);
    return data[index];
  }

  const double &operator()(int i, int j) const {
    std::size_t index = (static_cast<std::size_t>(i) * static_cast<std::size_t>(cols)) + static_cast<std::size_t>(j);
    return data[index];
  }
};

}  // namespace yakimov_i_mult_of_dense_matrices_fox_algorithm_seq
