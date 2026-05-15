#include "viderman_a_sparse_matrix_mult_crs_complex/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <thread>
#include <vector>

#include "util/include/util.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/common/include/common.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {

VidermanASparseMatrixMultCRSComplexSTL::VidermanASparseMatrixMultCRSComplexSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CRSMatrix();
}

void VidermanASparseMatrixMultCRSComplexSTL::ProcessRows(int start_row, int end_row,
                                                         std::vector<std::vector<Complex>> &local_values,
                                                         std::vector<std::vector<int>> &local_cols) const {
  const int cols_b = b_->cols;
  std::vector<Complex> accumulator(cols_b, Complex(0.0, 0.0));
  std::vector<int> marker(cols_b, -1);

  for (int i = start_row; i < end_row; ++i) {
    std::vector<int> current_row_indices;
    for (int j = a_->row_ptr[i]; j < a_->row_ptr[i + 1]; ++j) {
      int col_a = a_->col_indices[j];
      Complex val_a = a_->values[j];

      for (int k = b_->row_ptr[col_a]; k < b_->row_ptr[col_a + 1]; ++k) {
        int col_b = b_->col_indices[k];
        accumulator[col_b] += val_a * b_->values[k];
        if (marker[col_b] != i) {
          current_row_indices.push_back(col_b);
          marker[col_b] = i;
        }
      }
    }

    std::ranges::sort(current_row_indices);

    for (int idx : current_row_indices) {
      if (std::abs(accumulator[idx]) > kEpsilon) {
        local_values[i].push_back(accumulator[idx]);
        local_cols[i].push_back(idx);
      }
      accumulator[idx] = Complex(0.0, 0.0);
    }
  }
}

bool VidermanASparseMatrixMultCRSComplexSTL::ValidationImpl() {
  const auto &input = GetInput();
  const auto &a = std::get<0>(input);
  const auto &b = std::get<1>(input);
  return a.IsValid() && b.IsValid() && (a.cols == b.rows);
}

bool VidermanASparseMatrixMultCRSComplexSTL::PreProcessingImpl() {
  const auto &input = GetInput();
  a_ = &std::get<0>(input);
  b_ = &std::get<1>(input);
  return true;
}

bool VidermanASparseMatrixMultCRSComplexSTL::RunImpl() {
  const int rows_a = a_->rows;
  const int cols_b = b_->cols;
  const int num_threads = ppc::util::GetNumThreads();

  std::vector<std::jthread> threads;
  std::vector<std::vector<Complex>> local_values(rows_a);
  std::vector<std::vector<int>> local_cols(rows_a);

  int chunk = rows_a / num_threads;
  int remainder = rows_a % num_threads;
  int current_start = 0;

  for (int i = 0; i < num_threads; ++i) {
    int current_end = current_start + chunk + (i < remainder ? 1 : 0);
    if (current_start < current_end) {
      threads.emplace_back(&VidermanASparseMatrixMultCRSComplexSTL::ProcessRows, this, current_start, current_end,
                           std::ref(local_values), std::ref(local_cols));
    }
    current_start = current_end;
  }

  threads.clear();

  CRSMatrix &c = GetOutput();
  c.rows = rows_a;
  c.cols = cols_b;
  c.row_ptr.assign(rows_a + 1, 0);

  c.values.clear();
  c.col_indices.clear();

  for (int i = 0; i < rows_a; ++i) {
    c.row_ptr[i + 1] = c.row_ptr[i] + static_cast<int>(local_values[i].size());
    std::move(local_values[i].begin(), local_values[i].end(), std::back_inserter(c.values));
    std::move(local_cols[i].begin(), local_cols[i].end(), std::back_inserter(c.col_indices));
  }

  return true;
}

bool VidermanASparseMatrixMultCRSComplexSTL::PostProcessingImpl() {
  return GetOutput().IsValid();
}

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
