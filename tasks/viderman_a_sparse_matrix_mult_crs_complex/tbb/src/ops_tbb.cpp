#include "viderman_a_sparse_matrix_mult_crs_complex/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "oneapi/tbb/enumerable_thread_specific.h"
#include "oneapi/tbb/parallel_for.h"
#include "viderman_a_sparse_matrix_mult_crs_complex/common/include/common.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {
namespace {

void ProcessRow(const CRSMatrix &a, const CRSMatrix &b, int row, std::vector<Complex> &accumulator,
                std::vector<int> &marker, std::vector<int> &current_row_indices, std::vector<int> &dst_cols,
                std::vector<Complex> &dst_vals) {
  current_row_indices.clear();

  for (int j = a.row_ptr[row]; j < a.row_ptr[row + 1]; ++j) {
    const int col_a = a.col_indices[j];
    const Complex val_a = a.values[j];

    for (int k = b.row_ptr[col_a]; k < b.row_ptr[col_a + 1]; ++k) {
      const int col_b = b.col_indices[k];
      accumulator[col_b] += val_a * b.values[k];

      if (marker[col_b] != row) {
        current_row_indices.push_back(col_b);
        marker[col_b] = row;
      }
    }
  }

  std::ranges::sort(current_row_indices);

  dst_cols.clear();
  dst_vals.clear();
  dst_cols.reserve(current_row_indices.size());
  dst_vals.reserve(current_row_indices.size());

  for (const int idx : current_row_indices) {
    if (std::abs(accumulator[idx]) > kEpsilon) {
      dst_cols.push_back(idx);
      dst_vals.push_back(accumulator[idx]);
    }
    accumulator[idx] = Complex(0.0, 0.0);
  }
}

struct ThreadBuffers {
  std::vector<Complex> accumulator;
  std::vector<int> marker;
  std::vector<int> current_row_indices;
};

}  // namespace

VidermanASparseMatrixMultCRSComplexTBB::VidermanASparseMatrixMultCRSComplexTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CRSMatrix();
}

bool VidermanASparseMatrixMultCRSComplexTBB::ValidationImpl() {
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

bool VidermanASparseMatrixMultCRSComplexTBB::PreProcessingImpl() {
  const auto &input = GetInput();

  a_ = &std::get<0>(input);
  b_ = &std::get<1>(input);

  return true;
}

void VidermanASparseMatrixMultCRSComplexTBB::Multiply(const CRSMatrix &a, const CRSMatrix &b, CRSMatrix &c) {
  c.rows = a.rows;
  c.cols = b.cols;
  c.row_ptr.assign(a.rows + 1, 0);
  c.col_indices.clear();
  c.values.clear();

  std::vector<std::vector<int>> row_cols(a.rows);
  std::vector<std::vector<Complex>> row_vals(a.rows);

  oneapi::tbb::enumerable_thread_specific<ThreadBuffers> tls([&b] {
    ThreadBuffers buffers;
    buffers.accumulator.assign(b.cols, Complex(0.0, 0.0));
    buffers.marker.assign(b.cols, -1);
    return buffers;
  });

  oneapi::tbb::parallel_for(0, a.rows, [&](int i) {
    auto &buffers = tls.local();
    ProcessRow(a, b, i, buffers.accumulator, buffers.marker, buffers.current_row_indices, row_cols[i], row_vals[i]);
  });

  for (int i = 0; i < a.rows; ++i) {
    c.row_ptr[i + 1] = c.row_ptr[i] + static_cast<int>(row_cols[i].size());
  }

  c.col_indices.reserve(static_cast<std::size_t>(c.row_ptr[a.rows]));
  c.values.reserve(static_cast<std::size_t>(c.row_ptr[a.rows]));
  for (int i = 0; i < a.rows; ++i) {
    c.col_indices.insert(c.col_indices.end(), row_cols[i].begin(), row_cols[i].end());
    c.values.insert(c.values.end(), row_vals[i].begin(), row_vals[i].end());
  }
}

bool VidermanASparseMatrixMultCRSComplexTBB::RunImpl() {
  if (a_ == nullptr || b_ == nullptr) {
    return false;
  }

  CRSMatrix &c = GetOutput();
  Multiply(*a_, *b_, c);

  return true;
}

bool VidermanASparseMatrixMultCRSComplexTBB::PostProcessingImpl() {
  CRSMatrix &c = GetOutput();
  return c.IsValid();
}

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
