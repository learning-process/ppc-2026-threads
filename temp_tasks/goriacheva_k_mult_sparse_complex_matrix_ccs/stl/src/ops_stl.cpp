#include "goriacheva_k_mult_sparse_complex_matrix_ccs/stl/include/ops_stl.hpp"

#include <algorithm>
#include <ranges>
#include <thread>
#include <vector>

#include "goriacheva_k_mult_sparse_complex_matrix_ccs/common/include/common.hpp"

namespace goriacheva_k_mult_sparse_complex_matrix_ccs {

GoriachevaKMultSparseComplexMatrixCcsSTL::GoriachevaKMultSparseComplexMatrixCcsSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool GoriachevaKMultSparseComplexMatrixCcsSTL::ValidationImpl() {
  auto &[a, b] = GetInput();

  if (a.cols != b.rows) {
    return false;
  }

  if (a.col_ptr.empty() || b.col_ptr.empty()) {
    return false;
  }

  return true;
}

bool GoriachevaKMultSparseComplexMatrixCcsSTL::PreProcessingImpl() {
  GetOutput() = {};
  return true;
}

namespace {
void ProcessColumn(int j, const SparseMatrixCCS &a, const SparseMatrixCCS &b, std::vector<Complex> &values,
                   std::vector<int> &rows) {
  std::vector<Complex> accumulator(a.rows);
  std::vector<int> marker(a.rows, -1);
  std::vector<int> used_rows;

  for (int bi = b.col_ptr[j]; bi < b.col_ptr[j + 1]; bi++) {
    int k = b.row_ind[bi];
    Complex b_val = b.values[bi];

    for (int ai = a.col_ptr[k]; ai < a.col_ptr[k + 1]; ai++) {
      int i = a.row_ind[ai];

      if (marker[i] != j) {
        marker[i] = j;
        accumulator[i] = Complex(0.0, 0.0);
        used_rows.push_back(i);
      }

      accumulator[i] += a.values[ai] * b_val;
    }
  }

  std::ranges::sort(used_rows);

  for (int r : used_rows) {
    if (accumulator[r] != Complex(0.0, 0.0)) {
      rows.push_back(r);
      values.push_back(accumulator[r]);
    }
  }
}
}  // namespace

bool GoriachevaKMultSparseComplexMatrixCcsSTL::RunImpl() {
  auto &a = std::get<0>(GetInput());
  auto &b = std::get<1>(GetInput());
  auto &c = GetOutput();

  c.rows = a.rows;
  c.cols = b.cols;
  c.col_ptr.resize(c.cols + 1);

  std::vector<std::vector<Complex>> local_values(c.cols);
  std::vector<std::vector<int>> local_rows(c.cols);

  int num_threads = static_cast<int>(std::thread::hardware_concurrency());
  if (num_threads == 0) {
    num_threads = 1;
  }

  std::vector<std::thread> threads;

  auto worker = [&](int start, int end) {
    for (int j = start; j < end; j++) {
      ProcessColumn(j, a, b, local_values[j], local_rows[j]);
    }
  };

  int cols_per_thread = (b.cols + num_threads - 1) / num_threads;

  for (int thread_idx = 0; thread_idx < num_threads; thread_idx++) {
    int start = thread_idx * cols_per_thread;
    int end = std::min(start + cols_per_thread, b.cols);

    if (start >= end) {
      break;
    }

    threads.emplace_back(worker, start, end);
  }

  for (auto &th : threads) {
    th.join();
  }

  int nnz = 0;
  for (int j = 0; j < c.cols; j++) {
    c.col_ptr[j] = nnz;
    nnz += static_cast<int>(local_values[j].size());
  }
  c.col_ptr[c.cols] = nnz;

  c.values.reserve(nnz);
  c.row_ind.reserve(nnz);

  for (int j = 0; j < c.cols; j++) {
    c.values.insert(c.values.end(), local_values[j].begin(), local_values[j].end());
    c.row_ind.insert(c.row_ind.end(), local_rows[j].begin(), local_rows[j].end());
  }

  return true;
}

bool GoriachevaKMultSparseComplexMatrixCcsSTL::PostProcessingImpl() {
  return true;
}

}  // namespace goriacheva_k_mult_sparse_complex_matrix_ccs
