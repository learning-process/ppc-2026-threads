#include "ashihmin_d_mult_matr_crs/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <future>
#include <map>
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

bool AshihminDMultMatrCrsSTL::RunImpl() {
  const auto &matrix_a = GetInput().first;
  const auto &matrix_b = GetInput().second;
  auto &matrix_c = GetOutput();
  int rows_a = matrix_a.rows;

  std::vector<std::vector<int>> local_cols(rows_a);
  std::vector<std::vector<double>> local_vals(rows_a);

  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  std::vector<std::future<void>> futures;
  int chunk_size = (rows_a + num_threads - 1) / num_threads;

  for (unsigned int t = 0; t < num_threads; ++t) {
    int start_row = t * chunk_size;
    int end_row = std::min(start_row + chunk_size, rows_a);

    if (start_row >= end_row) {
      break;
    }

    futures.push_back(
        std::async(std::launch::async, [start_row, end_row, &matrix_a, &matrix_b, &local_cols, &local_vals]() {
      for (int i = start_row; i < end_row; ++i) {
        std::map<int, double> row_accumulator;

        for (int j = matrix_a.row_ptr[i]; j < matrix_a.row_ptr[i + 1]; ++j) {
          int col_a = matrix_a.col_index[j];
          double val_a = matrix_a.values[j];

          for (int k = matrix_b.row_ptr[col_a]; k < matrix_b.row_ptr[col_a + 1]; ++k) {
            int col_b = matrix_b.col_index[k];
            double val_b = matrix_b.values[k];
            row_accumulator[col_b] += val_a * val_b;
          }
        }

        for (auto it = row_accumulator.begin(); it != row_accumulator.end(); ++it) {
          if (std::abs(it->second) > 1e-15) {
            local_cols[i].push_back(it->first);
            local_vals[i].push_back(it->second);
          }
        }
      }
    }));
  }

  for (auto &f : futures) {
    f.get();
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
