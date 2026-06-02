#include "volkov_a_sparse_mat_mul_ccs/stl/include/ops_stl.hpp"

#include <cmath>
#include <cstddef>
#include <thread>
#include <tuple>
#include <vector>

#include "volkov_a_sparse_mat_mul_ccs/common/include/common.hpp"

namespace volkov_a_sparse_mat_mul_ccs {

namespace {

template <typename Func>
void RunParallelChunks(int total_size, const Func &worker_func) {
  if (total_size <= 0) {
    return;
  }

  int num_threads = static_cast<int>(std::thread::hardware_concurrency());
  if (num_threads == 0) {
    num_threads = 4;
  }

  int chunk_size = total_size / num_threads;
  int remainder = total_size % num_threads;

  std::vector<std::thread> threads;
  threads.reserve(static_cast<size_t>(num_threads));

  int current_start = 0;
  for (int i = 0; i < num_threads; ++i) {
    int current_end = current_start + chunk_size + (i < remainder ? 1 : 0);
    if (current_start < current_end) {
      threads.emplace_back(worker_func, current_start, current_end);
    }
    current_start = current_end;
  }

  for (auto &t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }
}

template <typename MatrixType>
void ProcessColumn(int col_idx, const MatrixType &matrix_a, const MatrixType &matrix_b,
                   std::vector<double> &col_accumulator, std::vector<int> &local_row_indices,
                   std::vector<double> &local_values) {
  int b_start = matrix_b.col_ptrs[col_idx];
  int b_end = matrix_b.col_ptrs[col_idx + 1];

  for (int k = b_start; k < b_end; ++k) {
    int b_row = matrix_b.row_indices[k];
    double b_val = matrix_b.values[k];

    int a_start = matrix_a.col_ptrs[b_row];
    int a_end = matrix_a.col_ptrs[b_row + 1];

    for (int idx = a_start; idx < a_end; ++idx) {
      int a_row = matrix_a.row_indices[idx];
      double a_val = matrix_a.values[idx];
      col_accumulator[a_row] += a_val * b_val;
    }
  }

  for (int i = 0; i < matrix_a.rows_count; ++i) {
    if (std::abs(col_accumulator[i]) > 1e-10) {
      local_row_indices.push_back(i);
      local_values.push_back(col_accumulator[i]);
    }
    col_accumulator[i] = 0.0;
  }
}

}  // namespace

VolkovASparseMatMulCcsStl::VolkovASparseMatMulCcsStl(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool VolkovASparseMatMulCcsStl::ValidationImpl() {
  const auto &matrix_a = std::get<0>(GetInput());
  const auto &matrix_b = std::get<1>(GetInput());
  return (matrix_a.cols_count == matrix_b.rows_count);
}

bool VolkovASparseMatMulCcsStl::PreProcessingImpl() {
  return true;
}

bool VolkovASparseMatMulCcsStl::RunImpl() {
  const auto &matrix_a = std::get<0>(GetInput());
  const auto &matrix_b = std::get<1>(GetInput());
  auto &matrix_c = GetOutput();

  matrix_c.rows_count = matrix_a.rows_count;
  matrix_c.cols_count = matrix_b.cols_count;

  std::vector<std::vector<int>> local_row_indices(matrix_b.cols_count);
  std::vector<std::vector<double>> local_values(matrix_b.cols_count);

  RunParallelChunks(matrix_b.cols_count, [&](int start, int end) {
    std::vector<double> col_accumulator(matrix_a.rows_count, 0.0);
    for (int j = start; j < end; ++j) {
      ProcessColumn(j, matrix_a, matrix_b, col_accumulator, local_row_indices[j], local_values[j]);
    }
  });

  matrix_c.col_ptrs.assign(matrix_c.cols_count + 1, 0);
  for (int j = 0; j < matrix_b.cols_count; ++j) {
    matrix_c.col_ptrs[j + 1] = matrix_c.col_ptrs[j] + static_cast<int>(local_row_indices[j].size());
  }

  matrix_c.non_zeros = matrix_c.col_ptrs.back();
  matrix_c.row_indices.resize(matrix_c.non_zeros);
  matrix_c.values.resize(matrix_c.non_zeros);

  RunParallelChunks(matrix_b.cols_count, [&](int start, int end) {
    for (int j = start; j < end; ++j) {
      int offset = matrix_c.col_ptrs[j];
      int current_col_size = static_cast<int>(local_row_indices[j].size());

      for (int k = 0; k < current_col_size; ++k) {
        matrix_c.row_indices[offset + k] = local_row_indices[j][k];
        matrix_c.values[offset + k] = local_values[j][k];
      }
    }
  });

  return true;
}

bool VolkovASparseMatMulCcsStl::PostProcessingImpl() {
  return true;
}

}  // namespace volkov_a_sparse_mat_mul_ccs
