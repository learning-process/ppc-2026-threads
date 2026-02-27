#include "potashnik_m_matrix_mult_complex/seq/include/ops_seq.hpp"

#include <chrono>
#include <thread>

#include "potashnik_m_matrix_mult_complex/common/include/common.hpp"

namespace potashnik_m_matrix_mult_complex {

PotashnikMMatrixMultComplexSEQ::PotashnikMMatrixMultComplexSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool PotashnikMMatrixMultComplexSEQ::ValidationImpl() {
  const auto &matrix_left = std::get<0>(GetInput());
  const auto &matrix_right = std::get<1>(GetInput());
  return matrix_left.width == matrix_right.height;
}

bool PotashnikMMatrixMultComplexSEQ::PreProcessingImpl() {
  return true;
}

bool PotashnikMMatrixMultComplexSEQ::RunImpl() {
  const auto &matrix_left = std::get<0>(GetInput());
  const auto &matrix_right = std::get<1>(GetInput());

  std::vector<Complex> val_left = matrix_left.val;
  std::vector<size_t> row_ind_left = matrix_left.row_ind;
  std::vector<size_t> col_ptr_left = matrix_left.col_ptr;
  size_t height_left = matrix_left.height;
  // size_t width_left = matrix_left.width;

  std::vector<Complex> val_right = matrix_right.val;
  std::vector<size_t> row_ind_right = matrix_right.row_ind;
  std::vector<size_t> col_ptr_right = matrix_right.col_ptr;
  // size_t height_right = matrix_right.height;
  size_t width_right = matrix_right.width;

  std::map<std::pair<size_t, size_t>, Complex> buffer;  

  for (size_t i = 0; i < matrix_left.Count(); i++) {
    size_t row_left = row_ind_left[i];
    size_t col_left = col_ptr_left[i];
    Complex left_val = val_left[i];

    for (size_t j = 0; j < matrix_right.Count(); j++) {
      size_t row_right = row_ind_right[j];
      size_t col_right = col_ptr_right[j];
      Complex right_val = val_right[j];

      if (col_left == row_right) {
        buffer[{row_left, col_right}] += left_val * right_val;
      }
    }
  }

  CCSMatrix matrix_res;

  matrix_res.width = width_right;
  matrix_res.height = height_left;
  for (const auto &[key, value] : buffer) {
    matrix_res.val.push_back(value);
    matrix_res.row_ind.push_back(key.first);
    matrix_res.col_ptr.push_back(key.second);
  }

  GetOutput() = matrix_res;
  return true;
}

bool PotashnikMMatrixMultComplexSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace potashnik_m_matrix_mult_complex
