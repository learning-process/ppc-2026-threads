#include "kapanova_s_sparse_matrix_mult_ccs/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <vector>

#include "kapanova_s_sparse_matrix_mult_ccs/common/include/common.hpp"

namespace kapanova_s_sparse_matrix_mult_ccs {

KapanovaSSparseMatrixMultCCSOMP::KapanovaSSparseMatrixMultCCSOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KapanovaSSparseMatrixMultCCSOMP::ValidationImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());

  if (a.cols != b.rows) {
    return false;
  }
  if (a.rows == 0 || a.cols == 0 || b.rows == 0 || b.cols == 0) {
    return false;
  }
  if (a.col_ptrs.size() != static_cast<size_t>(a.cols + 1)) {
    return false;
  }
  if (b.col_ptrs.size() != static_cast<size_t>(b.cols + 1)) {
    return false;
  }
  if (a.col_ptrs[0] != 0 || b.col_ptrs[0] != 0) {
    return false;
  }
  if (a.col_ptrs[a.cols] != a.nnz) {
    return false;
  }
  if (b.col_ptrs[b.cols] != b.nnz) {
    return false;
  }
  if (a.values.size() != static_cast<size_t>(a.nnz) || a.row_indices.size() != static_cast<size_t>(a.nnz)) {
    return false;
  }
  if (b.values.size() != static_cast<size_t>(b.nnz) || b.row_indices.size() != static_cast<size_t>(b.nnz)) {
    return false;
  }

  return true;
}

bool KapanovaSSparseMatrixMultCCSOMP::PreProcessingImpl() {
  return true;
}

void KapanovaSSparseMatrixMultCCSOMP::ProcessColumn(int j, const CCSMatrix &a, const CCSMatrix &b,
                                                    std::vector<std::vector<double>> &thread_accum,
                                                    std::vector<std::vector<bool>> &thread_row_mask,
                                                    std::vector<std::vector<size_t>> &thread_active_rows,
                                                    std::vector<std::vector<std::vector<size_t>>> &thread_col_rows,
                                                    std::vector<std::vector<std::vector<double>>> &thread_col_vals) {
  int tid = omp_get_thread_num();

  for (size_t k = b.col_ptrs[j]; k < b.col_ptrs[j + 1]; ++k) {
    size_t row_b = b.row_indices[k];
    double val_b = b.values[k];

    for (size_t zc = a.col_ptrs[row_b]; zc < a.col_ptrs[row_b + 1]; ++zc) {
      size_t i = a.row_indices[zc];
      double val_a = a.values[zc];

      thread_accum[tid][i] += val_a * val_b;
      if (!thread_row_mask[tid][i]) {
        thread_row_mask[tid][i] = true;
        thread_active_rows[tid].push_back(i);
      }
    }
  }

  std::sort(thread_active_rows[tid].begin(), thread_active_rows[tid].end());

  for (size_t i : thread_active_rows[tid]) {
    if (thread_accum[tid][i] != 0.0) {
      thread_col_rows[tid][j].push_back(i);
      thread_col_vals[tid][j].push_back(thread_accum[tid][i]);
    }
    thread_accum[tid][i] = 0.0;
    thread_row_mask[tid][i] = false;
  }
  thread_active_rows[tid].clear();
}

void KapanovaSSparseMatrixMultCCSOMP::ComputeColumnSizes(
    int num_threads, size_t cols, const std::vector<std::vector<std::vector<size_t>>> &thread_col_rows,
    std::vector<size_t> &col_sizes) {
  for (int tid = 0; tid < num_threads; ++tid) {
    for (size_t j = 0; j < cols; ++j) {
      col_sizes[j] += thread_col_rows[tid][j].size();
    }
  }
}

void KapanovaSSparseMatrixMultCCSOMP::MergeThreadResults(
    int num_threads, size_t cols, CCSMatrix &c, const std::vector<std::vector<std::vector<size_t>>> &thread_col_rows,
    const std::vector<std::vector<std::vector<double>>> &thread_col_vals) {
  std::vector<size_t> current_pos(cols, 0);

  for (int tid = 0; tid < num_threads; ++tid) {
    for (size_t j = 0; j < cols; ++j) {
      size_t start = c.col_ptrs[j] + current_pos[j];
      size_t num_elements = thread_col_rows[tid][j].size();

      for (size_t idx = 0; idx < num_elements; ++idx) {
        size_t pos = start + idx;
        c.row_indices[pos] = thread_col_rows[tid][j][idx];
        c.values[pos] = thread_col_vals[tid][j][idx];
      }
      current_pos[j] += num_elements;
    }
  }
}

bool KapanovaSSparseMatrixMultCCSOMP::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  OutType &c = GetOutput();

  c.rows = a.rows;
  c.cols = b.cols;
  c.col_ptrs.assign(c.cols + 1, 0);
  c.nnz = 0;

  int num_threads = omp_get_max_threads();

  std::vector<std::vector<double>> thread_accum(num_threads);
  std::vector<std::vector<bool>> thread_row_mask(num_threads);
  std::vector<std::vector<size_t>> thread_active_rows(num_threads);
  std::vector<std::vector<std::vector<size_t>>> thread_col_rows(num_threads);
  std::vector<std::vector<std::vector<double>>> thread_col_vals(num_threads);

#pragma omp parallel default(none) \
    shared(a, c, num_threads, thread_accum, thread_row_mask, thread_col_rows, thread_col_vals)
  {
    int tid = omp_get_thread_num();
    thread_accum[tid].assign(a.rows, 0.0);
    thread_row_mask[tid].assign(a.rows, false);
    thread_col_rows[tid].resize(c.cols);
    thread_col_vals[tid].resize(c.cols);
  }

#pragma omp parallel for schedule(dynamic) default(none) \
    shared(a, b, c, thread_accum, thread_row_mask, thread_active_rows, thread_col_rows, thread_col_vals)
  for (size_t j = 0; j < c.cols; ++j) {
    ProcessColumn(static_cast<int>(j), a, b, thread_accum, thread_row_mask, thread_active_rows, thread_col_rows,
                  thread_col_vals);
  }

  std::vector<size_t> col_sizes(c.cols, 0);
  ComputeColumnSizes(num_threads, c.cols, thread_col_rows, col_sizes);

  size_t offset = 0;
  for (size_t j = 0; j < c.cols; ++j) {
    c.col_ptrs[j] = offset;
    offset += col_sizes[j];
  }
  c.col_ptrs[c.cols] = offset;
  c.nnz = offset;

  c.values.resize(c.nnz);
  c.row_indices.resize(c.nnz);

  MergeThreadResults(num_threads, c.cols, c, thread_col_rows, thread_col_vals);

  return true;
}

bool KapanovaSSparseMatrixMultCCSOMP::PostProcessingImpl() {
  return true;
}

}  // namespace kapanova_s_sparse_matrix_mult_ccs
