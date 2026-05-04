// #include "potashnik_m_matrix_mult_complex/seq/include/ops_seq.hpp"
#include "potashnik_m_matrix_mult_complex/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include "potashnik_m_matrix_mult_complex/common/include/common.hpp"

namespace potashnik_m_matrix_mult_complex {

PotashnikMMatrixMultComplexOMP::PotashnikMMatrixMultComplexOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool PotashnikMMatrixMultComplexOMP::ValidationImpl() {
  const auto &matrix_left = std::get<0>(GetInput());
  const auto &matrix_right = std::get<1>(GetInput());
  return matrix_left.width == matrix_right.height;
}

bool PotashnikMMatrixMultComplexOMP::PreProcessingImpl() {
  return true;
}

bool PotashnikMMatrixMultComplexOMP::RunImpl() {
  const auto &matrix_left = std::get<0>(GetInput());
  const auto &matrix_right = std::get<1>(GetInput());

  const auto &val_left = matrix_left.val;
  const auto &row_ind_left = matrix_left.row_ind;
  const auto &col_ptr_left = matrix_left.col_ptr;
  size_t height_left = matrix_left.height;

  const auto &val_right = matrix_right.val;
  const auto &row_ind_right = matrix_right.row_ind;
  const auto &col_ptr_right = matrix_right.col_ptr;
  size_t width_right = matrix_right.width;

  int threads = omp_get_max_threads();
  std::vector<std::map<std::pair<size_t, size_t>, Complex>> local_buffers(threads);

#pragma omp parallel num_threads(threads) default(none)                                                              \
    shared(local_buffers, matrix_left, matrix_right, val_left, row_ind_left, col_ptr_left, val_right, row_ind_right, \
               col_ptr_right, width_right, height_left, threads)
  {
    int tid = omp_get_thread_num();
    auto &local_buffer = local_buffers[tid];
    size_t left_count = matrix_left.Count();
#pragma omp for schedule(static)
    for (size_t i = 0; i < left_count; i++) {
      size_t row_left = row_ind_left[i];
      size_t col_left = col_ptr_left[i];
      Complex left_val = val_left[i];

      for (size_t j = 0; j < matrix_right.Count(); j++) {
        size_t row_right = row_ind_right[j];
        size_t col_right = col_ptr_right[j];
        Complex right_val = val_right[j];

        if (col_left == row_right) {
          local_buffer[{row_left, col_right}] += left_val * right_val;
        }
      }
    }
  }

  std::map<std::pair<size_t, size_t>, Complex> buffer;
  for (const auto &local : local_buffers) {
    for (const auto &[key, value] : local) {
      buffer[key] += value;
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

bool PotashnikMMatrixMultComplexOMP::PostProcessingImpl() {
  return true;
}
}  // namespace potashnik_m_matrix_mult_complex
