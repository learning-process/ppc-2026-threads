#include "sinev_a_mult_matrix_fox_algorithm/seq/include/ops_seq.hpp"

#include <cstddef>
#include <vector>

#include "sinev_a_mult_matrix_fox_algorithm/common/include/common.hpp"

namespace sinev_a_mult_matrix_fox_algorithm_seq {

SinevAMultMatrixFoxAlgorithmSEQ::SinevAMultMatrixFoxAlgorithmSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool SinevAMultMatrixFoxAlgorithmSEQ::ValidationImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();

  return matrix_size > 0 && matrix_a.size() == matrix_size * matrix_size &&
         matrix_b.size() == matrix_size * matrix_size;
}

bool SinevAMultMatrixFoxAlgorithmSEQ::PreProcessingImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  GetOutput() = std::vector<double>(matrix_size * matrix_size, 0.0);
  return true;
}

bool SinevAMultMatrixFoxAlgorithmSEQ::RunImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();

  size_t n = matrix_size;
  auto &output = GetOutput();

  if (output.size() != n * n) {
    output.assign(n * n, 0.0);
  }

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < n; ++k) {
        sum += matrix_a[(i * n) + k] * matrix_b[(k * n) + j];
      }
      output[(i * n) + j] = sum;
    }
  }

  return true;
}

bool SinevAMultMatrixFoxAlgorithmSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace sinev_a_mult_matrix_fox_algorithm_seq
