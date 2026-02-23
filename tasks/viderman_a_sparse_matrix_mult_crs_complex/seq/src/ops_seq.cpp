#include "viderman_a_sparse_matrix_mult_crs_complex/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace viderman_a_sparse_matrix_mult_crs_complex {

VidermanASparseMatrixMultCRSComplexSEQ::VidermanASparseMatrixMultCRSComplexSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CRSMatrix();
  A_ = nullptr;
  B_ = nullptr;
}

bool VidermanASparseMatrixMultCRSComplexSEQ::ValidationImpl() {
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

bool VidermanASparseMatrixMultCRSComplexSEQ::PreProcessingImpl() {
  const auto &input = GetInput();

  A_ = &std::get<0>(input);
  B_ = &std::get<1>(input);

  return true;
}

void VidermanASparseMatrixMultCRSComplexSEQ::Multiply(const CRSMatrix &A, const CRSMatrix &B, CRSMatrix &C) {
  C.rows = A.rows;
  C.cols = B.cols;
  C.row_ptr.assign(A.rows + 1, 0);
  C.col_indices.clear();
  C.values.clear();

  std::vector<Complex> accumulator(B.cols, Complex(0.0, 0.0));
  std::vector<int> marker(B.cols, -1);
  std::vector<int> current_row_indices;

  for (int i = 0; i < A.rows; ++i) {
    current_row_indices.clear();

    for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) {
      int col_A = A.col_indices[j];
      Complex val_A = A.values[j];

      for (int k = B.row_ptr[col_A]; k < B.row_ptr[col_A + 1]; ++k) {
        int col_B = B.col_indices[k];
        accumulator[col_B] += val_A * B.values[k];

        if (marker[col_B] != i) {
          current_row_indices.push_back(col_B);
          marker[col_B] = i;
        }
      }
    }

    std::sort(current_row_indices.begin(), current_row_indices.end());

    C.row_ptr[i + 1] = C.row_ptr[i];
    for (int idx : current_row_indices) {
      if (std::abs(accumulator[idx]) > EPSILON) {
        C.values.push_back(accumulator[idx]);
        C.col_indices.push_back(idx);
        ++C.row_ptr[i + 1];
      }
      accumulator[idx] = Complex(0.0, 0.0);
    }
  }
}

bool VidermanASparseMatrixMultCRSComplexSEQ::RunImpl() {
  if (!A_ || !B_) {
    return false;
  }

  CRSMatrix &C = GetOutput();
  Multiply(*A_, *B_, C);

  return true;
}

bool VidermanASparseMatrixMultCRSComplexSEQ::PostProcessingImpl() {
  CRSMatrix &C = GetOutput();
  return C.IsValid();
}

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
