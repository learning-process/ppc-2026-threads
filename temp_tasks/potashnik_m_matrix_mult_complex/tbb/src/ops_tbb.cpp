#include "potashnik_m_matrix_mult_complex/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include "potashnik_m_matrix_mult_complex/common/include/common.hpp"

namespace potashnik_m_matrix_mult_complex {

PotashnikMMatrixMultComplexTBB::PotashnikMMatrixMultComplexTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool PotashnikMMatrixMultComplexTBB::ValidationImpl() {
  const auto &matrix_left = std::get<0>(GetInput());
  const auto &matrix_right = std::get<1>(GetInput());
  return matrix_left.width == matrix_right.height;
}

bool PotashnikMMatrixMultComplexTBB::PreProcessingImpl() {
  return true;
}

bool PotashnikMMatrixMultComplexTBB::RunImpl() {
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

  size_t left_count = matrix_left.Count();
  size_t right_count = matrix_right.Count();

  using Key = std::pair<size_t, size_t>;
  using LocalMap = std::map<Key, Complex>;

  tbb::enumerable_thread_specific<LocalMap> local_buffers;

  const size_t grain_i = 64;
  const size_t grain_j = 256;

  tbb::parallel_for(tbb::blocked_range2d<size_t>(0, left_count, grain_i, 0, right_count, grain_j),
                    [&](const tbb::blocked_range2d<size_t> &r) {
    auto &local_buffer = local_buffers.local();

    for (size_t i = r.rows().begin(); i != r.rows().end(); ++i) {
      size_t row_left = row_ind_left[i];
      size_t col_left = col_ptr_left[i];
      Complex left_val = val_left[i];

      for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
        if (col_left == row_ind_right[j]) {
          local_buffer[{row_left, col_ptr_right[j]}] += left_val * val_right[j];
        }
      }
    }
  });

  std::map<Key, Complex> buffer;
  for (const auto &local : local_buffers) {
    for (const auto &[key, value] : local) {
      buffer[key] += value;
    }
  }

  CCSMatrix matrix_res;
  matrix_res.width = width_right;
  matrix_res.height = height_left;

  matrix_res.val.reserve(buffer.size());
  matrix_res.row_ind.reserve(buffer.size());
  matrix_res.col_ptr.reserve(buffer.size());

  for (const auto &[key, value] : buffer) {
    matrix_res.val.push_back(value);
    matrix_res.row_ind.push_back(key.first);
    matrix_res.col_ptr.push_back(key.second);
  }

  GetOutput() = matrix_res;
  return true;
}

bool PotashnikMMatrixMultComplexTBB::PostProcessingImpl() {
  return true;
}

}  // namespace potashnik_m_matrix_mult_complex
