#include "luzan_e_double_sparse_matrix_mult/stl/include/ops_stl.hpp"

#include "luzan_e_double_sparse_matrix_mult/common/include/common.hpp"
// #include "util/include/util.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

namespace luzan_e_double_sparse_matrix_mult {
void LuzanEDoubleSparseMatrixMultSTL::AccumulateColumn(const SparseMatrix &a, const SparseMatrix &b, unsigned b_col,
                                                       std::vector<double> &tmp_col) {
  unsigned b_rows_start = b.col_index[b_col];
  unsigned b_rows_end = b.col_index[b_col + 1];

  for (unsigned b_pos = b_rows_start; b_pos < b_rows_end; b_pos++) {
    double b_val = b.value[b_pos];
    unsigned b_row = b.row[b_pos];
    unsigned a_rows_start = a.col_index[b_row];
    unsigned a_rows_end = a.col_index[b_row + 1];

    for (unsigned a_pos = a_rows_start; a_pos < a_rows_end; a_pos++) {
      tmp_col[a.row[a_pos]] += a.value[a_pos] * b_val;
    }
  }
}

void LuzanEDoubleSparseMatrixMultSTL::CollectNonZeros(const std::vector<double> &tmp_col, unsigned b_col,
                                                      std::vector<std::vector<double>> &values_per_col,
                                                      std::vector<std::vector<unsigned>> &rows_per_col) {
  for (unsigned i = 0; i < static_cast<unsigned>(tmp_col.size()); i++) {
    if (fabs(tmp_col[i]) > kEPS) {
      values_per_col[b_col].push_back(tmp_col[i]);
      rows_per_col[b_col].push_back(i);
    }
  }
}

void LuzanEDoubleSparseMatrixMultSTL::ProcessColumn(const SparseMatrix &a, const SparseMatrix &b, unsigned b_col,
                                                    std::vector<std::vector<double>> &values_per_col,
                                                    std::vector<std::vector<unsigned>> &rows_per_col) {
  std::vector<double> tmp_col(a.rows, 0.0);
  AccumulateColumn(a, b, b_col, tmp_col);
  CollectNonZeros(tmp_col, b_col, values_per_col, rows_per_col);
}

void LuzanEDoubleSparseMatrixMultSTL::AssembleResult(SparseMatrix &c, unsigned cols,
                                                     const std::vector<std::vector<double>> &values_per_col,
                                                     const std::vector<std::vector<unsigned>> &rows_per_col) {
  c.col_index.push_back(0);
  for (unsigned j = 0; j < cols; j++) {
    for (size_t k = 0; k < values_per_col[j].size(); k++) {
      c.value.push_back(values_per_col[j][k]);
      c.row.push_back(rows_per_col[j][k]);
    }
    c.col_index.push_back(static_cast<unsigned>(c.value.size()));
  }
}

SparseMatrix LuzanEDoubleSparseMatrixMultSTL::CalcProdSTL(const SparseMatrix &a, const SparseMatrix &b) {
  SparseMatrix c(a.rows, b.cols);

  std::vector<std::vector<double>> values_per_col(b.cols);
  std::vector<std::vector<unsigned>> rows_per_col(b.cols);

  const unsigned num_threads = std::thread::hardware_concurrency();
  const unsigned chunk = (b.cols + num_threads - 1) / num_threads;

  auto worker = [&](unsigned thread_id) {
    unsigned start = thread_id * chunk;
    unsigned end = std::min(start + chunk, static_cast<unsigned>(b.cols));
    for (unsigned b_col = start; b_col < end; b_col++) {
      ProcessColumn(a, b, b_col, values_per_col, rows_per_col);
    }
  };

  std::vector<std::thread> threads(num_threads);
  for (unsigned thrd = 0; thrd < num_threads; thrd++) {
    threads[thrd] = std::thread(worker, thrd);
  }
  for (auto &thrd : threads) {
    thrd.join();
  }

  AssembleResult(c, b.cols, values_per_col, rows_per_col);
  return c;
}

LuzanEDoubleSparseMatrixMultSTL::LuzanEDoubleSparseMatrixMultSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  // GetOutput() = 0;
}

bool LuzanEDoubleSparseMatrixMultSTL::ValidationImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  return a.GetCols() == b.GetRows() && a.GetCols() != 0 && a.GetRows() != 0 && b.GetCols() != 0;
}

bool LuzanEDoubleSparseMatrixMultSTL::PreProcessingImpl() {
  return true;
}

bool LuzanEDoubleSparseMatrixMultSTL::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());

  GetOutput() = CalcProdSTL(a, b);
  return true;
}

bool LuzanEDoubleSparseMatrixMultSTL::PostProcessingImpl() {
  return true;
}

}  // namespace luzan_e_double_sparse_matrix_mult
