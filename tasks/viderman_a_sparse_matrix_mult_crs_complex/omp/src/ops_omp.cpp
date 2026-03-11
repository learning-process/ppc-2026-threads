#include "viderman_a_sparse_matrix_mult_crs_complex/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "viderman_a_sparse_matrix_mult_crs_complex/common/include/common.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {

VidermanASparseMatrixMultCRSComplexOMP::VidermanASparseMatrixMultCRSComplexOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CRSMatrix();
}

bool VidermanASparseMatrixMultCRSComplexOMP::ValidationImpl() {
  const auto &input = GetInput();
  const auto &A = std::get<0>(input);
  const auto &B = std::get<1>(input);

  if (!A.IsValid() || !B.IsValid()) {
    return false;
  }

  if (A.cols != B.rows) {
    return false;
  }

  return true;
}

bool VidermanASparseMatrixMultCRSComplexOMP::PreProcessingImpl() {
  const auto &input = GetInput();

  A_ = &std::get<0>(input);
  B_ = &std::get<1>(input);

  return true;
}

void VidermanASparseMatrixMultCRSComplexOMP::Multiply(const CRSMatrix &A, const CRSMatrix &B, CRSMatrix &C) {
  C.rows = A.rows;
  C.cols = B.cols;
  C.row_ptr.assign(A.rows + 1, 0);
  C.col_indices.clear();
  C.values.clear();

  std::vector<std::vector<int>> row_cols(A.rows);
  std::vector<std::vector<Complex>> row_vals(A.rows);

#pragma omp parallel shared(A, B, row_cols, row_vals)
  {
    std::vector<Complex> accumulator(B.cols, Complex(0.0, 0.0));
    std::vector<int> marker(B.cols, -1);
    std::vector<int> current_row_indices;

#pragma omp for schedule(static)
    for (int i = 0; i < A.rows; ++i) {
      current_row_indices.clear();

      for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) {
        const int col_a = A.col_indices[j];
        const Complex val_a = A.values[j];

        for (int k = B.row_ptr[col_a]; k < B.row_ptr[col_a + 1]; ++k) {
          const int col_b = B.col_indices[k];
          accumulator[col_b] += val_a * B.values[k];

          if (marker[col_b] != i) {
            current_row_indices.push_back(col_b);
            marker[col_b] = i;
          }
        }
      }

      std::ranges::sort(current_row_indices);

      auto &dst_cols = row_cols[i];
      auto &dst_vals = row_vals[i];
      dst_cols.clear();
      dst_vals.clear();
      dst_cols.reserve(current_row_indices.size());
      dst_vals.reserve(current_row_indices.size());

      for (const int idx : current_row_indices) {
        if (std::abs(accumulator[idx]) > kEpsilon) {
          dst_cols.push_back(idx);
          dst_vals.push_back(accumulator[idx]);
        }
        accumulator[idx] = Complex(0.0, 0.0);
      }
    }
  }

  for (int i = 0; i < A.rows; ++i) {
    C.row_ptr[i + 1] = C.row_ptr[i] + static_cast<int>(row_cols[i].size());
  }

  C.col_indices.reserve(static_cast<std::size_t>(C.row_ptr[A.rows]));
  C.values.reserve(static_cast<std::size_t>(C.row_ptr[A.rows]));
  for (int i = 0; i < A.rows; ++i) {
    C.col_indices.insert(C.col_indices.end(), row_cols[i].begin(), row_cols[i].end());
    C.values.insert(C.values.end(), row_vals[i].begin(), row_vals[i].end());
  }
}

bool VidermanASparseMatrixMultCRSComplexOMP::RunImpl() {
  if (A_ == nullptr || B_ == nullptr) {
    return false;
  }

  CRSMatrix &C = GetOutput();
  Multiply(*A_, *B_, C);

  return true;
}

bool VidermanASparseMatrixMultCRSComplexOMP::PostProcessingImpl() {
  CRSMatrix &C = GetOutput();
  return C.IsValid();
}

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
