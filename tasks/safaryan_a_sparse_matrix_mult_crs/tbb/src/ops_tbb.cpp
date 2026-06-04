#include "safaryan_a_sparse_matrix_mult_crs/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/parallel_for.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "safaryan_a_sparse_matrix_mult_crs/common/include/common.hpp"
#include "util/include/util.hpp"

namespace safaryan_a_sparse_matrix_mult_crs {

SafaryanATaskTBB::SafaryanATaskTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = SparseMatrixCCS();
}

bool SafaryanATaskTBB::IsMatrixValid(const SparseMatrixCCS &matrix) {
  if (matrix.rows < 0 || matrix.cols < 0) {
    return false;
  }
  if (matrix.col_ptrs.size() != static_cast<size_t>(matrix.cols) + 1) {
    return false;
  }
  if (matrix.values.size() != matrix.row_indices.size()) {
    return false;
  }

  if (matrix.col_ptrs.empty() || matrix.col_ptrs[0] != 0) {
    return false;
  }

  const int total_elements = static_cast<int>(matrix.values.size());
  if (matrix.col_ptrs[matrix.cols] != total_elements) {
    return false;
  }

  for (size_t i = 0; i < matrix.col_ptrs.size() - 1; ++i) {
    if (matrix.col_ptrs[i] > matrix.col_ptrs[i + 1] || matrix.col_ptrs[i] < 0) {
      return false;
    }
  }

  for (size_t i = 0; i < matrix.row_indices.size(); ++i) {
    if (matrix.row_indices[i] < 0 || matrix.row_indices[i] >= matrix.rows) {
      return false;
    }
  }

  return true;
}

bool SafaryanATaskTBB::ValidationImpl() {
  const auto &[a, b] = GetInput();

  if (!IsMatrixValid(a) || !IsMatrixValid(b)) {
    return false;
  }
  if (a.cols != b.rows) {
    return false;
  }

  return true;
}

bool SafaryanATaskTBB::PreProcessingImpl() {
  const auto &[a, b] = GetInput();
  GetOutput() = SparseMatrixCCS(a.rows, b.cols);
  return true;
}

namespace {
void ProcessColumn(const SparseMatrixCCS &a, const SparseMatrixCCS &b, int j, std::vector<double> &temp) {
  for (int b_idx = b.col_ptrs[j]; b_idx < b.col_ptrs[j + 1]; ++b_idx) {
    const int k = b.row_indices[b_idx];
    const double b_val = b.values[b_idx];
    for (int a_idx = a.col_ptrs[k]; a_idx < a.col_ptrs[k + 1]; ++a_idx) {
      const int i = a.row_indices[a_idx];
      const double a_val = a.values[a_idx];
      temp[i] += a_val * b_val;
    }
  }
}

SparseMatrixCCS BuildColumnRange(const SparseMatrixCCS &a, const SparseMatrixCCS &b, int begin_col, int end_col) {
  SparseMatrixCCS part(a.rows, end_col - begin_col);
  std::vector<double> temp(a.rows, 0.0);

  for (int j = begin_col; j < end_col; ++j) {
    part.col_ptrs[j - begin_col] = static_cast<int>(part.values.size());
    ProcessColumn(a, b, j, temp);

    for (int i = 0; i < a.rows; ++i) {
      if (std::abs(temp[i]) > 1e-10) {
        part.values.push_back(temp[i]);
        part.row_indices.push_back(i);
        temp[i] = 0.0;
      }
    }
  }

  part.col_ptrs[end_col - begin_col] = static_cast<int>(part.values.size());
  return part;
}
}  // namespace

SparseMatrixCCS SafaryanATaskTBB::MultiplyMatrices(const SparseMatrixCCS &a, const SparseMatrixCCS &b) {
  SparseMatrixCCS result(a.rows, b.cols);
  if (b.cols == 0) {
    return result;
  }

  const int chunks_count = std::max(1, std::min(ppc::util::GetNumThreads(), b.cols));
  std::vector<SparseMatrixCCS> parts(chunks_count);
  std::vector<int> begin_cols(chunks_count + 1, 0);

  for (int i = 0; i <= chunks_count; ++i) {
    begin_cols[i] = (b.cols * i) / chunks_count;
  }

  tbb::parallel_for(0, chunks_count, [&](int chunk_id) {
    parts[chunk_id] = BuildColumnRange(a, b, begin_cols[chunk_id], begin_cols[chunk_id + 1]);
  });

  for (int part_id = 0; part_id < chunks_count; ++part_id) {
    const SparseMatrixCCS &part = parts[part_id];
    const int column_offset = begin_cols[part_id];
    for (int local_col = 0; local_col < part.cols; ++local_col) {
      result.col_ptrs[column_offset + local_col] = static_cast<int>(result.values.size());
      const int begin = part.col_ptrs[local_col];
      const int end = part.col_ptrs[local_col + 1];
      result.values.insert(result.values.end(), part.values.begin() + begin, part.values.begin() + end);
      result.row_indices.insert(result.row_indices.end(), part.row_indices.begin() + begin,
                                part.row_indices.begin() + end);
    }
  }

  result.col_ptrs[b.cols] = static_cast<int>(result.values.size());
  return result;
}

bool SafaryanATaskTBB::RunImpl() {
  const auto &[a, b] = GetInput();
  GetOutput() = SparseMatrixCCS(a.rows, b.cols);
  GetOutput() = MultiplyMatrices(a, b);
  return true;
}

bool SafaryanATaskTBB::PostProcessingImpl() {
  return true;
}

}  // namespace safaryan_a_sparse_matrix_mult_crs
