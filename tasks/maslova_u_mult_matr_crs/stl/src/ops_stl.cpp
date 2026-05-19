#include "maslova_u_mult_matr_crs/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

#include "maslova_u_mult_matr_crs/common/include/common.hpp"
#include "util/include/util.hpp"

namespace maslova_u_mult_matr_crs {

MaslovaUMultMatrSTL::MaslovaUMultMatrSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool MaslovaUMultMatrSTL::ValidationImpl() {
  const auto &input = GetInput();
  const auto &matrix_a = std::get<0>(input);
  const auto &matrix_b = std::get<1>(input);

  if (matrix_a.cols != matrix_b.rows || matrix_a.rows <= 0 || matrix_b.cols <= 0) {
    return false;
  }
  if (matrix_a.row_ptr.size() != static_cast<size_t>(matrix_a.rows) + 1) {
    return false;
  }
  if (matrix_b.row_ptr.size() != static_cast<size_t>(matrix_b.rows) + 1) {
    return false;
  }
  return true;
}

bool MaslovaUMultMatrSTL::PreProcessingImpl() {
  const auto &matrix_a = std::get<0>(GetInput());
  const auto &matrix_b = std::get<1>(GetInput());
  auto &matrix_c = GetOutput();
  matrix_c.rows = matrix_a.rows;
  matrix_c.cols = matrix_b.cols;
  return true;
}

int MaslovaUMultMatrSTL::GetRowNNZ(int i, const CRSMatrix &a, const CRSMatrix &b, std::vector<int> &marker) {
  int row_nnz = 0;
  for (int j = a.row_ptr[i]; j < a.row_ptr[i + 1]; ++j) {
    int col_a = a.col_ind[j];
    for (int k = b.row_ptr[col_a]; k < b.row_ptr[col_a + 1]; ++k) {
      int col_b = b.col_ind[k];
      if (marker[col_b] != i) {
        marker[col_b] = i;
        row_nnz++;
      }
    }
  }
  return row_nnz;
}

void MaslovaUMultMatrSTL::FillRowValues(int i, const CRSMatrix &a, const CRSMatrix &b, CRSMatrix &c,
                                        std::vector<double> &acc, std::vector<int> &marker, std::vector<int> &used) {
  used.clear();
  for (int j = a.row_ptr[i]; j < a.row_ptr[i + 1]; ++j) {
    int col_a = a.col_ind[j];
    double val_a = a.values[j];
    for (int k = b.row_ptr[col_a]; k < b.row_ptr[col_a + 1]; ++k) {
      int col_b = b.col_ind[k];
      if (marker[col_b] != i) {
        marker[col_b] = i;
        used.push_back(col_b);
        acc[col_b] = val_a * b.values[k];
      } else {
        acc[col_b] += val_a * b.values[k];
      }
    }
  }

  if (!used.empty()) {
    std::ranges::sort(used);
    int write_pos = c.row_ptr[i];
    for (int col : used) {
      c.values[write_pos] = acc[col];
      c.col_ind[write_pos] = col;
      write_pos++;
      acc[col] = 0.0;
    }
  }
}

void MaslovaUMultMatrSTL::ComputeRowPtrSTL(int rows_a, int cols_b, int num_threads) {
  const auto &matrix_a = std::get<0>(GetInput());
  const auto &matrix_b = std::get<1>(GetInput());
  auto &matrix_c = GetOutput();

  std::vector<std::thread> threads;
  int chunk = rows_a / num_threads;

  for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    int start = thread_idx * chunk;
    int end = (thread_idx == num_threads - 1) ? rows_a : (thread_idx + 1) * chunk;
    if (start < end) {
      threads.emplace_back([&, start, end, cols_b]() {
        std::vector<int> marker(cols_b, -1);
        for (int i = start; i < end; ++i) {
          matrix_c.row_ptr[i + 1] = GetRowNNZ(i, matrix_a, matrix_b, marker);
        }
      });
    }
  }
  for (auto &thread_item : threads) {
    thread_item.join();
  }
}

void MaslovaUMultMatrSTL::ComputeValuesSTL(int rows_a, int cols_b, int num_threads) {
  const auto &matrix_a = std::get<0>(GetInput());
  const auto &matrix_b = std::get<1>(GetInput());
  auto &matrix_c = GetOutput();

  std::vector<std::thread> threads;
  int chunk = rows_a / num_threads;

  for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    int start = thread_idx * chunk;
    int end = (thread_idx == num_threads - 1) ? rows_a : (thread_idx + 1) * chunk;
    if (start < end) {
      threads.emplace_back([&, start, end, cols_b]() {
        std::vector<double> acc(cols_b, 0.0);
        std::vector<int> marker(cols_b, -1);
        std::vector<int> used;
        used.reserve(cols_b);
        for (int i = start; i < end; ++i) {
          FillRowValues(i, matrix_a, matrix_b, matrix_c, acc, marker, used);
        }
      });
    }
  }
  for (auto &thread_item : threads) {
    thread_item.join();
  }
}

bool MaslovaUMultMatrSTL::RunImpl() {
  const auto &matrix_a = std::get<0>(GetInput());
  const auto &matrix_b = std::get<1>(GetInput());
  auto &matrix_c = GetOutput();

  int rows_a = matrix_a.rows;
  int cols_b = matrix_b.cols;
  matrix_c.row_ptr.assign(static_cast<size_t>(rows_a) + 1, 0);

  int num_threads = std::max(1, ppc::util::GetNumThreads());
  num_threads = std::min(num_threads, rows_a);

  ComputeRowPtrSTL(rows_a, cols_b, num_threads);

  for (int i = 0; i < rows_a; ++i) {
    matrix_c.row_ptr[static_cast<size_t>(i) + 1] += matrix_c.row_ptr[i];
  }

  matrix_c.values.resize(matrix_c.row_ptr[rows_a]);
  matrix_c.col_ind.resize(matrix_c.row_ptr[rows_a]);

  ComputeValuesSTL(rows_a, cols_b, num_threads);

  return true;
}

bool MaslovaUMultMatrSTL::PostProcessingImpl() {
  return true;
}

}  // namespace maslova_u_mult_matr_crs
