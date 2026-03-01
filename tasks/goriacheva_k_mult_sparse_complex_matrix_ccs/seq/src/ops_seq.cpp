#include "goriacheva_k_mult_sparse_complex_matrix_ccs/seq/include/ops_seq.hpp"

#include "goriacheva_k_mult_sparse_complex_matrix_ccs/common/include/common.hpp"
#include "util/include/util.hpp"

namespace goriacheva_k_mult_sparse_complex_matrix_ccs {

GoriachevaKMultSparseComplexMatrixCcsSEQ::GoriachevaKMultSparseComplexMatrixCcsSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool GoriachevaKMultSparseComplexMatrixCcsSEQ::ValidationImpl() {
  auto &[A, B] = GetInput();

  if (A.cols != B.rows) {
    return false;
  }

  if (A.col_ptr.empty() || B.col_ptr.empty()) {
    return false;
  }

  return true;
}

bool GoriachevaKMultSparseComplexMatrixCcsSEQ::PreProcessingImpl() {
  GetOutput() = {};
  return true;
}

bool GoriachevaKMultSparseComplexMatrixCcsSEQ::RunImpl() {
  auto &[A, B] = GetInput();
  auto &C = GetOutput();

  C.rows = A.rows;
  C.cols = B.cols;
  C.col_ptr.resize(C.cols + 1);

  std::vector<Complex> values;
  std::vector<int> rows;

  std::vector<Complex> accumulator(A.rows);
  std::vector<int> marker(A.rows, -1);
  std::vector<int> used_rows;

  for (int j = 0; j < B.cols; j++) {
    C.col_ptr[j] = values.size();
    used_rows.clear();

    for (int bi = B.col_ptr[j]; bi < B.col_ptr[j + 1]; bi++) {
      int k = B.row_ind[bi];
      Complex b_val = B.values[bi];

      for (int ai = A.col_ptr[k]; ai < A.col_ptr[k + 1]; ai++) {
        int i = A.row_ind[ai];

        if (marker[i] != j) {
          marker[i] = j;
          accumulator[i] = Complex(0.0, 0.0);
          used_rows.push_back(i);
        }

        accumulator[i] += A.values[ai] * b_val;
      }
    }

    std::sort(used_rows.begin(), used_rows.end());

    for (int r : used_rows) {
      if (accumulator[r] != Complex(0.0, 0.0)) {
        rows.push_back(r);
        values.push_back(accumulator[r]);
      }
    }
  }

  C.col_ptr[C.cols] = values.size();
  C.values = std::move(values);
  C.row_ind = std::move(rows);

  return true;
}

bool GoriachevaKMultSparseComplexMatrixCcsSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace goriacheva_k_mult_sparse_complex_matrix_ccs
