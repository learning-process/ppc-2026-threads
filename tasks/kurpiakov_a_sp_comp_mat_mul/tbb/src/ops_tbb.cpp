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

struct ThreadLocalRowBuffers {
  explicit ThreadLocalRowBuffers(int cols) : row_acc(cols), row_used(cols, 0) {}

  std::vector<ComplexD> row_acc;
  std::vector<char> row_used;
  std::vector<int> used_cols;
};

void ComputeRow(const SparseMatrix &a, const SparseMatrix &b, int row_idx, ThreadLocalRowBuffers &buffers,
                std::vector<ComplexD> &out_values, std::vector<int> &out_cols) {
  buffers.used_cols.clear();

  for (int ja = a.row_ptr[row_idx]; ja < a.row_ptr[row_idx + 1]; ++ja) {
    const int ka = a.col_indices[ja];
    const ComplexD &a_val = a.values[ja];

    for (int jb = b.row_ptr[ka]; jb < b.row_ptr[ka + 1]; ++jb) {
      const int cb = b.col_indices[jb];
      const ComplexD &b_val = b.values[jb];

      if (buffers.row_used[cb] == 0) {
        buffers.row_used[cb] = 1;
        buffers.row_acc[cb] = ComplexD();
        buffers.used_cols.push_back(cb);
      }

      buffers.row_acc[cb] += a_val * b_val;
    }
  }

  std::ranges::sort(buffers.used_cols);

  out_values.clear();
  out_cols.clear();
  out_values.reserve(buffers.used_cols.size());
  out_cols.reserve(buffers.used_cols.size());

  for (int c : buffers.used_cols) {
    out_values.push_back(buffers.row_acc[c]);
    out_cols.push_back(c);
    buffers.row_used[c] = 0;
  }
}

SparseMatrix BuildResultFromRows(int rows, int cols, const std::vector<std::vector<ComplexD>> &row_values,
                                 const std::vector<std::vector<int>> &row_cols) {
  SparseMatrix result(rows, cols);

  for (int i = 0; i < rows; ++i) {
    result.row_ptr[i + 1] = result.row_ptr[i] + static_cast<int>(row_values[i].size());
  }

  const auto total_nnz = static_cast<std::size_t>(result.row_ptr[rows]);

  result.values.resize(total_nnz);
  result.col_indices.resize(total_nnz);

  tbb::parallel_for(tbb::blocked_range<int>(0, rows), [&](const tbb::blocked_range<int> &range) {
    for (int i = range.begin(); i < range.end(); ++i) {
      const auto offset = static_cast<std::vector<ComplexD>::difference_type>(result.row_ptr[i]);
      std::copy(row_values[i].begin(), row_values[i].end(), result.values.begin() + offset);
      std::copy(row_cols[i].begin(), row_cols[i].end(), result.col_indices.begin() + offset);
    }
  });

  return result;
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

  tbb::enumerable_thread_specific<ThreadLocalRowBuffers> tls_buffers([&] { return ThreadLocalRowBuffers(cols); });

  tbb::parallel_for(tbb::blocked_range<int>(0, rows), [&](const tbb::blocked_range<int> &range) {
    auto &buffers = tls_buffers.local();
    for (int i = range.begin(); i < range.end(); ++i) {
      ComputeRow(a, b, i, buffers, row_values[i], row_cols[i]);
    }
  });

  GetOutput() = BuildResultFromRows(rows, cols, row_values, row_cols);
  return true;
}

bool KurpiakovACRSMatMulTBB::PostProcessingImpl() {
  return true;
}

}  // namespace kurpiakov_a_sp_comp_mat_mul
