#include "tsyplakov_k_mul_double_crs_matrix/stl/include/ops_stl.hpp"

#include <cmath>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include "tsyplakov_k_mul_double_crs_matrix/common/include/common.hpp"
#include "util/include/util.hpp"

namespace tsyplakov_k_mul_double_crs_matrix {

namespace {

void ComputeRow(const SparseMatrixCRS& a, const SparseMatrixCRS& b, int row, std::vector<double>& values,
                std::vector<int>& cols) {
  std::unordered_map<int, double> acc;

  for (int idx_a = a.row_ptr[row]; idx_a < a.row_ptr[row + 1]; ++idx_a) {
    const int k = a.col_index[idx_a];
    const double val_a = a.values[idx_a];

    for (int idx_b = b.row_ptr[k]; idx_b < b.row_ptr[k + 1]; ++idx_b) {
      const int j = b.col_index[idx_b];
      acc[j] += val_a * b.values[idx_b];
    }
  }

  values.reserve(acc.size());
  cols.reserve(acc.size());

  for (const auto& [col, val] : acc) {
    if (std::fabs(val) > 1e-12) {
      cols.push_back(col);
      values.push_back(val);
    }
  }
}

}  // namespace

TsyplakovKTestTaskSTL::TsyplakovKTestTaskSTL(const InType& in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TsyplakovKTestTaskSTL::ValidationImpl() {
  const auto& input = GetInput();
  return input.a.cols == input.b.rows;
}

bool TsyplakovKTestTaskSTL::PreProcessingImpl() {
  return true;
}

bool TsyplakovKTestTaskSTL::RunImpl() {
  const auto& input = GetInput();
  const auto& a = input.a;
  const auto& b = input.b;

  const int rows = a.rows;
  const int num_threads = ppc::util::GetNumThreads();

  std::vector<std::vector<double>> row_values(rows);
  std::vector<std::vector<int>> row_cols(rows);

  std::vector<std::thread> threads(num_threads);

  const int chunk = (rows + num_threads - 1) / num_threads;

  for (int t = 0; t < num_threads; ++t) {
    threads[t] = std::thread([&, t]() {
      const int start = t * chunk;
      const int end = std::min(start + chunk, rows);

      for (int i = start; i < end; ++i) {
        ComputeRow(a, b, i, row_values[i], row_cols[i]);
      }
    });
  }

  for (auto& th : threads) {
    th.join();
  }

  SparseMatrixCRS c(a.rows, b.cols);

  for (int i = 0; i < c.rows; ++i) {
    c.row_ptr[i + 1] = c.row_ptr[i] + static_cast<int>(row_values[i].size());
  }

  const int nnz = c.row_ptr[c.rows];

  c.values.reserve(nnz);
  c.col_index.reserve(nnz);

  for (int i = 0; i < c.rows; ++i) {
    c.values.insert(c.values.end(), row_values[i].begin(), row_values[i].end());

    c.col_index.insert(c.col_index.end(), row_cols[i].begin(), row_cols[i].end());
  }

  GetOutput() = std::move(c);

  return true;
}

bool TsyplakovKTestTaskSTL::PostProcessingImpl() {
  return true;
}

}  // namespace tsyplakov_k_mul_double_crs_matrix
