#pragma once

#include <cstddef>
#include <vector>

#include "task/include/task.hpp"

namespace akhmetov_daniil_strassen_dense_double_seq {

using Matrix = std::vector<double>;
using InType = Matrix;
using OutType = Matrix;
using TestType = int;
using BaseTask = ppc::task::Task<InType, OutType>;

namespace format {

inline size_t GetN(const InType &input) {
  return static_cast<size_t>(input[0]);
}

inline size_t MatrixSize(const InType &input) {
  return GetN(input) * GetN(input);
}

inline Matrix GetA(const InType &input) {
  size_t n = GetN(input);
  return {input.begin() + 1, input.begin() + 1 + static_cast<std::ptrdiff_t>(n * n)};
}

inline Matrix GetB(const InType &input) {
  size_t n = GetN(input);
  return {input.begin() + 1 + static_cast<std::ptrdiff_t>(n * n),
          input.begin() + 1 + static_cast<std::ptrdiff_t>(2 * n * n)};
}

}  // namespace format
}  // namespace akhmetov_daniil_strassen_dense_double_seq
