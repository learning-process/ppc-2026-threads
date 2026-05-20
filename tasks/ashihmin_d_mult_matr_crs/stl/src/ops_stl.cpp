#include "ashihmin_d_mult_matr_crs/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <future>
#include <map>
#include <thread>
#include <vector>

#include "ashihmin_d_mult_matr_crs/common/include/common.hpp"

namespace ashihmin_d_mult_matr_crs {

AshihminDMultMatrCrsSTL::AshihminDMultMatrCrsSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool AshihminDMultMatrCrsSTL::ValidationImpl() {
  return GetInput().first.cols == GetInput().second.rows;
}

bool AshihminDMultMatrCrsSTL::PreProcessingImpl() {
  auto &matrix_c = GetOutput();
  matrix_c.rows = GetInput().first.rows;
  matrix_c.cols = GetInput().second.cols;
  matrix_c.row_ptr.assign(matrix_c.rows + 1, 0);
  matrix_c.values.clear();
  matrix_c.col_index.clear();
  return true;
}

void AshihminDMultMatrCrsSTL::MultiplyRow(int row_idx, const CRSMatrix &matrix_a, const CRSMatrix &matrix_b,
                                          std::vector<int> &row_cols, std::vector<double> &row_vals) {
  std::map<int, double> row_accumulator;
  for (int j = matrix_a.row_ptr[row_idx]; j < matrix_a.row_ptr[row_idx + 1]; ++j) {
    int col_a = matrix_a.col_index[j];
    double val_a = matrix_a.values[j];
    for (int k = matrix_b.row_ptr[col_a]; k < matrix_b.row_ptr[col_a + 1]; ++k) {
      row_accumulator[matrix_b.col_index[k]] += val_a * matrix_b.values[k];
    }
  }
  for (const auto &entry : row_accumulator) {
    if (std::abs(entry.second) > 1e-15) {
      row_cols.push_back(entry.first);
      row_vals.push_back(entry.second);
    }
  }
}

bool AshihminDMultMatrCrsSTL::RunImpl() {
  const auto &matrix_a = GetInput().first;
  const auto &matrix_b = GetInput().second;
  auto &matrix_c = GetOutput();
  int rows_a = matrix_a.rows;

  std::vector<std::vector<int>> local_cols(rows_a);
  std::vector<std::vector<double>> local_vals(rows_a);

  unsigned int hardware_threads = std::thread::hardware_concurrency();
  int num_threads = (hardware_threads == 0) ? 2 : static_cast<int>(hardware_threads);

  int chunk_size = (rows_a + num_threads - 1) / num_threads;
  std::vector<std::future<void>> futures;

  for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    int start_row = thread_idx * chunk_size;
    int end_row = std::min(start_row + chunk_size, rows_a);

    if (start_row >= end_row) {
      break;
    }

    futures.push_back(
        std::async(std::launch::async, [start_row, end_row, &matrix_a, &matrix_b, &local_cols, &local_vals] {
      for (int i = start_row; i < end_row; ++i) {
        MultiplyRow(i, matrix_a, matrix_b, local_cols[i], local_vals[i]);
      }
    }));
  }

  for (auto &fut : futures) {
    fut.get();
  }

  for (int i = 0; i < rows_a; ++i) {
    matrix_c.col_index.insert(matrix_c.col_index.end(), local_cols[i].begin(), local_cols[i].end());
    matrix_c.values.insert(matrix_c.values.end(), local_vals[i].begin(), local_vals[i].end());
    matrix_c.row_ptr[i + 1] = static_cast<int>(matrix_c.values.size());
  }

  return true;
}

bool AshihminDMultMatrCrsSTL::PostProcessingImpl() {
  return true;
}

}  // namespace ashihmin_d_mult_matr_crs
