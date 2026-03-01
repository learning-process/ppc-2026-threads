#include "goriacheva_k_mult_sparse_complex_matrix_ccs/seq/include/ops_seq.hpp"

#include <algorithm>
#include <utility>
#include <vector>

#include "goriacheva_k_mult_sparse_complex_matrix_ccs/common/include/common.hpp"

namespace goriacheva_k_mult_sparse_complex_matrix_ccs {

GoriachevaKMultSparseComplexMatrixCcsSEQ::GoriachevaKMultSparseComplexMatrixCcsSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool GoriachevaKMultSparseComplexMatrixCcsSEQ::ValidationImpl() {
  auto &[a, b] = GetInput();

  if (a.cols != b.rows) {
    return false;
  }

  if (a.col_ptr.empty() || b.col_ptr.empty()) {
    return false;
  }

  return true;
}

bool GoriachevaKMultSparseComplexMatrixCcsSEQ::PreProcessingImpl() {
  GetOutput() = {};
  return true;
}

bool GoriachevaKMultSparseComplexMatrixCcsSEQ::RunImpl() {
  auto &[a, b] = GetInput();
  auto &c = GetOutput();

  c.rows = a.rows;
  c.cols = b.cols;
  c.col_ptr.resize(c.cols + 1);

  std::vector<Complex> values;
  std::vector<int> rows;

  std::vector<Complex> accumulator(a.rows);
  std::vector<int> marker(a.rows, -1);
  std::vector<int> used_rows;

  for (int j = 0; j < b.cols; j++) {
    c.col_ptr[j] = static_cast<int>(values.size());
    used_rows.clear();

    for (int bi = b.col_ptr[j]; bi < b.col_ptr[j + 1]; bi++) {
      int k = b.row_ind[bi];
      Complex b_val = b.values[bi];

      for (int ai = a.col_ptr[k]; ai < a.col_ptr[k + 1]; ai++) {
        int i = a.row_ind[ai];

        if (marker[i] != j) {
          marker[i] = j;
          accumulator[i] = Complex(0.0, 0.0);
          used_rows.push_back(i);
        }

        accumulator[i] += a.values[ai] * b_val;
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

  c.col_ptr[c.cols] = static_cast<int>(values.size());
  c.values = std::move(values);
  c.row_ind = std::move(rows);

  return true;
}

bool GoriachevaKMultSparseComplexMatrixCcsSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace goriacheva_k_mult_sparse_complex_matrix_ccs
