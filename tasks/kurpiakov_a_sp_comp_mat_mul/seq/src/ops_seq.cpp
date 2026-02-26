#include "kurpiakov_a_sp_comp_mat_mul/seq/include/ops_seq.hpp"

#include "kurpiakov_a_sp_comp_mat_mul/common/include/common.hpp"

namespace kurpiakov_a_sp_comp_mat_mul {

namespace {

bool ValidateCSR(const SparseMatrix &m) {
  if (m.rows <= 0 || m.cols <= 0) {
    return false;
  }
  if (static_cast<int>(m.row_ptr.size()) != m.rows + 1) {
    return false;
  }
  if (m.row_ptr[0] != 0) {
    return false;
  }
  if (std::cmp_not_equal(m.values.size(), m.row_ptr[m.rows])) {
    return false;
  }
  if (m.col_indices.size() != m.values.size()) {
    return false;
  }
  for (int i = 0; i < m.rows; ++i) {
    for (int j = m.row_ptr[i]; j < m.row_ptr[i + 1]; ++j) {
      if (m.col_indices[j] < 0 || m.col_indices[j] >= m.cols) {
        return false;
      }
    }
  }
  return true;
}

}  // namespace

KurpiskovACRSMatMulSEQ::KurpiskovACRSMatMulSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = SparseMatrix();
}

bool KurpiskovACRSMatMulSEQ::ValidationImpl() {
  const auto &[a, b] = GetInput();

  if (!ValidateCSR(a) || !ValidateCSR(b)) {
    return false;
  }

  return a.cols == b.rows;
}

bool KurpiskovACRSMatMulSEQ::PreProcessingImpl() {
  return true;
}

bool KurpiskovACRSMatMulSEQ::RunImpl() {
  const auto &[a, b] = GetInput();
  GetOutput() = a.Multiply(b);
  return true;
}

bool KurpiskovACRSMatMulSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace kurpiakov_a_sp_comp_mat_mul
