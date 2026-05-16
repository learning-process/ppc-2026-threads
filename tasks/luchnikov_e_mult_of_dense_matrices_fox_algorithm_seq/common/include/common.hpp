#pragma once
#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace luchnikov_e_mult_of_dense_matrices_fox_algorithm_seq {

using InType = int;
using OutType = double;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

struct DenseMatrix {
  int rows = 0;
  int cols = 0;
  std::vector<double> values{};

  double &At(int r, int c) {
    return values[static_cast<std::size_t>(r) * cols + c];
  }

  const double &At(int r, int c) const {
    return values[static_cast<std::size_t>(r) * cols + c];
  }

  bool IsEmpty() const {
    return values.empty();
  }
  bool IsSquare() const {
    return rows == cols;
  }
};

}  // namespace luchnikov_e_mult_of_dense_matrices_fox_algorithm_seq
