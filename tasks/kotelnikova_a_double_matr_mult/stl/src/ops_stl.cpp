#include "kotelnikova_a_double_matr_mult/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

#include "kotelnikova_a_double_matr_mult/common/include/common.hpp"

namespace kotelnikova_a_double_matr_mult {

KotelnikovaATaskSTL::KotelnikovaATaskSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = SparseMatrixCCS();
}

bool KotelnikovaATaskSTL::IsMatrixValid(const SparseMatrixCCS &matrix) {
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

bool KotelnikovaATaskSTL::ValidationImpl() {
  const auto &[a, b] = GetInput();

  if (!IsMatrixValid(a) || !IsMatrixValid(b)) {
    return false;
  }
  if (a.cols != b.rows) {
    return false;
  }

  return true;
}

bool KotelnikovaATaskSTL::PreProcessingImpl() {
  const auto &[a, b] = GetInput();
  GetOutput() = SparseMatrixCCS(a.rows, b.cols);
  return true;
}

namespace {

std::vector<double> ComputeColumn(const SparseMatrixCCS &a, const SparseMatrixCCS &b, int col_idx) {
  std::vector<double> temp(static_cast<size_t>(a.rows), 0.0);

  for (int b_idx = b.col_ptrs[col_idx]; b_idx < b.col_ptrs[col_idx + 1]; ++b_idx) {
    const int k = b.row_indices[b_idx];
    const double b_val = b.values[b_idx];

    for (int a_idx = a.col_ptrs[k]; a_idx < a.col_ptrs[k + 1]; ++a_idx) {
      const int i = a.row_indices[a_idx];
      temp[static_cast<size_t>(i)] += a.values[a_idx] * b_val;
    }
  }

  return temp;
}

int CountNonZero(const std::vector<double> &column, double epsilon) {
  int count = 0;
  for (double val : column) {
    if (std::abs(val) > epsilon) {
      ++count;
    }
  }
  return count;
}

void FillColumn(const std::vector<double> &column, double epsilon, std::vector<int> &row_indices,
                std::vector<double> &values, int start_pos) {
  int pos = start_pos;
  for (size_t i = 0; i < column.size(); ++i) {
    if (std::abs(column[i]) > epsilon) {
      row_indices[static_cast<size_t>(pos)] = static_cast<int>(i);
      values[static_cast<size_t>(pos)] = column[i];
      ++pos;
    }
  }
}

}  // namespace

SparseMatrixCCS KotelnikovaATaskSTL::MultiplyMatrices(const SparseMatrixCCS &a, const SparseMatrixCCS &b) {
  SparseMatrixCCS result(a.rows, b.cols);

  const double epsilon = 1e-10;
  const unsigned int hardware_threads = std::thread::hardware_concurrency();
  const int num_threads = (hardware_threads > 0) ? static_cast<int>(hardware_threads) : 1;
  const int cols_per_thread = (b.cols + num_threads - 1) / num_threads;

  std::vector<int> col_start(b.cols, 0);
  std::vector<std::thread> threads;

  for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    int start_col = thread_idx * cols_per_thread;
    int end_col = std::min(start_col + cols_per_thread, b.cols);

    threads.emplace_back([&, start_col, end_col]() {
      for (int j = start_col; j < end_col; ++j) {
        std::vector<double> column = ComputeColumn(a, b, j);
        col_start[j] = CountNonZero(column, epsilon);
      }
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }

  std::vector<int> col_ptr(b.cols + 1, 0);
  for (int j = 0; j < b.cols; ++j) {
    col_ptr[j + 1] = col_ptr[j] + col_start[j];
  }

  const int total_nnz = col_ptr[b.cols];
  result.values.resize(static_cast<size_t>(total_nnz));
  result.row_indices.resize(static_cast<size_t>(total_nnz));
  result.col_ptrs = col_ptr;

  threads.clear();
  for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    int start_col = thread_idx * cols_per_thread;
    int end_col = std::min(start_col + cols_per_thread, b.cols);

    threads.emplace_back([&, start_col, end_col]() {
      for (int j = start_col; j < end_col; ++j) {
        std::vector<double> column = ComputeColumn(a, b, j);
        FillColumn(column, epsilon, result.row_indices, result.values, col_ptr[j]);
      }
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }

  return result;
}

bool KotelnikovaATaskSTL::RunImpl() {
  const auto &[a, b] = GetInput();
  GetOutput() = MultiplyMatrices(a, b);
  return true;
}

bool KotelnikovaATaskSTL::PostProcessingImpl() {
  return true;
}

}  // namespace kotelnikova_a_double_matr_mult
