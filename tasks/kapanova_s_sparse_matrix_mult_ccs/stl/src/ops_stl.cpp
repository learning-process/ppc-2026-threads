#include "kapanova_s_sparse_matrix_mult_ccs/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <thread>
#include <vector>

#include "kapanova_s_sparse_matrix_mult_ccs/common/include/common.hpp"

namespace kapanova_s_sparse_matrix_mult_ccs {

KapanovaSSparseMatrixMultCCSSTL::KapanovaSSparseMatrixMultCCSSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KapanovaSSparseMatrixMultCCSSTL::ValidationImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());

  if (a.cols != b.rows) {
    return false;
  }
  if (a.rows == 0 || a.cols == 0 || b.rows == 0 || b.cols == 0) {
    return false;
  }
  if (a.col_ptrs.size() != static_cast<size_t>(a.cols + 1)) {
    return false;
  }
  if (b.col_ptrs.size() != static_cast<size_t>(b.cols + 1)) {
    return false;
  }
  if (a.col_ptrs[0] != 0 || b.col_ptrs[0] != 0) {
    return false;
  }
  if (a.col_ptrs[a.cols] != a.nnz) {
    return false;
  }
  if (b.col_ptrs[b.cols] != b.nnz) {
    return false;
  }
  if (a.values.size() != static_cast<size_t>(a.nnz) || a.row_indices.size() != static_cast<size_t>(a.nnz)) {
    return false;
  }
  if (b.values.size() != static_cast<size_t>(b.nnz) || b.row_indices.size() != static_cast<size_t>(b.nnz)) {
    return false;
  }

  return true;
}

bool KapanovaSSparseMatrixMultCCSSTL::PreProcessingImpl() {
  return true;
}

bool KapanovaSSparseMatrixMultCCSSTL::PostProcessingImpl() {
  return true;
}

namespace {

struct ThreadLocalData {
  std::vector<double> accum;
  std::vector<char> mask;
  std::vector<size_t> active;

  explicit ThreadLocalData(size_t rows) {
    accum.assign(rows, 0.0);
    mask.assign(rows, 0);
    active.reserve(rows / 20);
  }

  void ClearActive() {
    for (size_t i : active) {
      accum[i] = 0.0;
      mask[i] = 0;
    }
    active.clear();
  }
};

void ProcessColumn(int j, const CCSMatrix &a, const CCSMatrix &b, std::vector<size_t> &out_rows,
                   std::vector<double> &out_vals, ThreadLocalData &local) {
  local.ClearActive();

  for (size_t k = b.col_ptrs[j]; k < b.col_ptrs[j + 1]; ++k) {
    size_t row_b = b.row_indices[k];
    double val_b = b.values[k];

    for (size_t zc = a.col_ptrs[row_b]; zc < a.col_ptrs[row_b + 1]; ++zc) {
      size_t i = a.row_indices[zc];
      double val_a = a.values[zc];

      local.accum[i] += val_a * val_b;
      if (local.mask[i] == 0) {
        local.mask[i] = 1;
        local.active.push_back(i);
      }
    }
  }

  std::ranges::sort(local.active);

  for (size_t i : local.active) {
    if (local.accum[i] != 0.0) {
      out_rows.push_back(i);
      out_vals.push_back(local.accum[i]);
    }
  }
}

void WorkerTask(const CCSMatrix &a, const CCSMatrix &b, int start, int end, std::vector<std::vector<size_t>> &rows,
                std::vector<std::vector<double>> &vals) {
  ThreadLocalData local(a.rows);

  for (int j = start; j < end; ++j) {
    ProcessColumn(j, a, b, rows[j], vals[j], local);
  }
}

}  // namespace

bool KapanovaSSparseMatrixMultCCSSTL::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  OutType &c = GetOutput();

  c.rows = a.rows;
  c.cols = b.cols;
  c.col_ptrs.assign(c.cols + 1, 0);
  c.nnz = 0;

  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  std::vector<std::vector<size_t>> temp_rows(c.cols);
  std::vector<std::vector<double>> temp_vals(c.cols);
  std::vector<std::thread> threads;

  size_t chunk = (c.cols + num_threads - 1) / num_threads;

  for (unsigned int thread_id = 0; thread_id < num_threads; ++thread_id) {
    size_t start = thread_id * chunk;
    size_t end = std::min(start + chunk, c.cols);
    if (start >= c.cols) {
      break;
    }

    threads.emplace_back(WorkerTask, std::cref(a), std::cref(b), static_cast<int>(start), static_cast<int>(end),
                         std::ref(temp_rows), std::ref(temp_vals));
  }

  for (auto &th : threads) {
    if (th.joinable()) {
      th.join();
    }
  }

  size_t offset = 0;
  for (size_t j = 0; j < c.cols; ++j) {
    c.col_ptrs[j] = offset;
    offset += temp_rows[j].size();
  }
  c.col_ptrs[c.cols] = offset;
  c.nnz = offset;

  c.values.resize(c.nnz);
  c.row_indices.resize(c.nnz);

  for (size_t j = 0; j < c.cols; ++j) {
    size_t start = c.col_ptrs[j];
    size_t n = temp_rows[j].size();
    for (size_t idx = 0; idx < n; ++idx) {
      c.row_indices[start + idx] = temp_rows[j][idx];
      c.values[start + idx] = temp_vals[j][idx];
    }
  }

  return true;
}

}  // namespace kapanova_s_sparse_matrix_mult_ccs
