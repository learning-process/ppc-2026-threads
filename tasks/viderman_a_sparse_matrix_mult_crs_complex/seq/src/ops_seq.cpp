#include "viderman_a_sparse_matrix_mult_crs_complex/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "viderman_a_sparse_matrix_mult_crs_complex/common/include/common.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {

VidermanASparseMatrixMultCRSComplexSEQ::VidermanASparseMatrixMultCRSComplexSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CRSMatrix();
}

bool VidermanASparseMatrixMultCRSComplexSEQ::ValidationImpl() {
  const auto &input = GetInput();
  const auto &a = std::get<0>(input);
  const auto &b = std::get<1>(input);

  if (!a.IsValid() || !b.IsValid()) {
    return false;
  }

  if (a.cols != b.rows) {
    return false;
  }

  return true;
}

bool VidermanASparseMatrixMultCRSComplexSEQ::PreProcessingImpl() {
  const auto &input = GetInput();

  a_ = &std::get<0>(input);
  b_ = &std::get<1>(input);

  return true;
}

void VidermanASparseMatrixMultCRSComplexSEQ::Multiply(const CRSMatrix &a, const CRSMatrix &b, CRSMatrix &c) {
  c.rows = a.rows;
  c.cols = b.cols;
  c.row_ptr.assign(a.rows + 1, 0);
  c.col_indices.clear();
  c.values.clear();

  std::vector<Complex> accumulator(b.cols, Complex(0.0, 0.0));
  std::vector<int> marker(b.cols, -1);
  std::vector<int> current_row_indices;

  for (int i = 0; i < a.rows; ++i) {
    current_row_indices.clear();

    for (int j = a.row_ptr[i]; j < a.row_ptr[i + 1]; ++j) {
      int col_a = a.col_indices[j];
      Complex val_a = a.values[j];

      for (int k = b.row_ptr[col_a]; k < b.row_ptr[col_a + 1]; ++k) {
        int col_b = b.col_indices[k];
        accumulator[col_b] += val_a * b.values[k];

        if (marker[col_b] != i) {
          current_row_indices.push_back(col_b);
          marker[col_b] = i;
        }
      }
    }

    std::ranges::sort(current_row_indices);

    c.row_ptr[i + 1] = c.row_ptr[i];
    for (int idx : current_row_indices) {
      if (std::abs(accumulator[idx]) > kEpsilon) {
        c.values.push_back(accumulator[idx]);
        c.col_indices.push_back(idx);
        ++c.row_ptr[i + 1];
      }
      accumulator[idx] = Complex(0.0, 0.0);
    }
  }
}

bool VidermanASparseMatrixMultCRSComplexSEQ::RunImpl() {
  if (a_ == nullptr || b_ == nullptr) {
    return false;
  }

  CRSMatrix &c = GetOutput();
  Multiply(*a_, *b_, c);

  return true;
}

bool VidermanASparseMatrixMultCRSComplexSEQ::PostProcessingImpl() {
  CRSMatrix &c = GetOutput();
  return c.IsValid();
}

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
