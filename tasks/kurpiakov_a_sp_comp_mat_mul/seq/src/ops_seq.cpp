#include "kurpiakov_a_sp_comp_mat_mul/seq/include/ops_seq.hpp"

namespace kurpiakov_a_sp_comp_mat_mul {

KurpiskovACRSMatMulSEQ::KurpiskovACRSMatMulSEQ(const InType& in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = SparseMatrix();
}

bool KurpiskovACRSMatMulSEQ::ValidationImpl() {
  const auto& [a, b] = GetInput();

  if (a._rows <= 0 || a._cols <= 0 || b._rows <= 0 || b._cols <= 0) {
    return false;
  }

  if (a._cols != b._rows) {
    return false;
  }

  for (int i = 0; i < a._rows; ++i) {
    for (int j = a._row_ptr[i]; j < a._row_ptr[i + 1]; ++j) {
      if (a._col_indices[j] < 0 || a._col_indices[j] >= a._cols) {
        return false;
      }
    }
  }

  if (static_cast<int>(b._row_ptr.size()) != b._rows + 1) {
    return false;
  }

  if (b._row_ptr[0] != 0) {
    return false;
  }

  for (int i = 0; i < b._rows; ++i) {
    for (int j = b._row_ptr[i]; j < b._row_ptr[i + 1]; ++j) {
      if (b._col_indices[j] < 0 || b._col_indices[j] >= b._cols) {
        return false;
      }
    }
  }

  return true;
}

bool KurpiskovACRSMatMulSEQ::PreProcessingImpl() {
  return true;
}

bool KurpiskovACRSMatMulSEQ::RunImpl() {
  const auto& [a, b] = GetInput();
  GetOutput() = a.Multiply(b);
  return true;
}

bool KurpiskovACRSMatMulSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace kurpiakov_a_sp_comp_mat_mul
