#include "kurpiakov_a_sp_comp_mat_mul/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

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

KurpiakovACRSMatMulTBB::KurpiakovACRSMatMulTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = SparseMatrix();
}

bool KurpiakovACRSMatMulTBB::ValidationImpl() {
  const auto &[a, b] = GetInput();

  if (!ValidateCSR(a) || !ValidateCSR(b)) {
    return false;
  }

  return a.cols == b.rows;
}

bool KurpiakovACRSMatMulTBB::PreProcessingImpl() {
  return true;
}

bool KurpiakovACRSMatMulTBB::RunImpl() {
  const auto &[a, b] = GetInput();
  const int rows = a.rows;
  const int cols = b.cols;

  std::vector<std::vector<ComplexD>> row_values(rows);
  std::vector<std::vector<int>> row_cols(rows);

  tbb::parallel_for(tbb::blocked_range<int>(0, rows), [&](const tbb::blocked_range<int> &range) {
    std::vector<ComplexD> row_acc(cols);
    std::vector<char> row_used(cols, 0);
    std::vector<int> used_cols;

    for (int i = range.begin(); i < range.end(); ++i) {
      used_cols.clear();

      for (int ja = a.row_ptr[i]; ja < a.row_ptr[i + 1]; ++ja) {
        const int ka = a.col_indices[ja];
        const ComplexD &a_val = a.values[ja];

        for (int jb = b.row_ptr[ka]; jb < b.row_ptr[ka + 1]; ++jb) {
          const int cb = b.col_indices[jb];
          const ComplexD &b_val = b.values[jb];

          if (row_used[cb] == 0) {
            row_used[cb] = 1;
            row_acc[cb] = ComplexD();
            used_cols.push_back(cb);
          }
          row_acc[cb] += a_val * b_val;
        }
      }

      std::ranges::sort(used_cols);

      auto &current_row_values = row_values[i];
      auto &current_row_cols = row_cols[i];
      current_row_values.clear();
      current_row_cols.clear();
      current_row_values.reserve(used_cols.size());
      current_row_cols.reserve(used_cols.size());

      for (int c : used_cols) {
        current_row_values.push_back(row_acc[c]);
        current_row_cols.push_back(c);
        row_used[c] = 0;
      }
    }
  });

  SparseMatrix result(rows, cols);
  std::size_t total_nnz = 0;
  for (int i = 0; i < rows; ++i) {
    total_nnz += row_values[i].size();
  }

  result.values.reserve(total_nnz);
  result.col_indices.reserve(total_nnz);

  for (int i = 0; i < rows; ++i) {
    result.values.insert(result.values.end(), row_values[i].begin(), row_values[i].end());
    result.col_indices.insert(result.col_indices.end(), row_cols[i].begin(), row_cols[i].end());
    result.row_ptr[i + 1] = static_cast<int>(result.values.size());
  }

  GetOutput() = std::move(result);
  return true;
}

bool KurpiakovACRSMatMulTBB::PostProcessingImpl() {
  return true;
}

}  // namespace kurpiakov_a_sp_comp_mat_mul
