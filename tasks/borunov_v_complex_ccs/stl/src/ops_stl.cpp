#include "borunov_v_complex_ccs/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include "borunov_v_complex_ccs/common/include/common.hpp"

namespace borunov_v_complex_ccs {

namespace {

void WorkerThread(int thread_id, int num_threads, int num_cols, const SparseMatrix &a, const SparseMatrix &b,
                  std::vector<std::complex<double>> &thread_val, std::vector<int> &thread_row_idx,
                  std::vector<int> &thread_col_ptr) {
  int start_col = (num_cols * thread_id) / num_threads;
  int end_col = (num_cols * (thread_id + 1)) / num_threads;

  int num_cols_thread = end_col - start_col;
  thread_col_ptr.assign(num_cols_thread + 1, 0);

  std::vector<std::complex<double>> col_accumulator(a.num_rows, {0.0, 0.0});
  std::vector<int> non_zero_indices;
  std::vector<bool> is_non_zero(a.num_rows, false);

  int current_nnz = 0;
  for (int j = start_col; j < end_col; ++j) {
    for (int b_idx = b.col_ptrs[j]; b_idx < b.col_ptrs[j + 1]; ++b_idx) {
      int p = b.row_indices[b_idx];
      std::complex<double> b_val = b.values[b_idx];

      for (int a_idx = a.col_ptrs[p]; a_idx < a.col_ptrs[p + 1]; ++a_idx) {
        int i = a.row_indices[a_idx];
        std::complex<double> a_val = a.values[a_idx];

        if (!is_non_zero[i]) {
          is_non_zero[i] = true;
          non_zero_indices.push_back(i);
        }
        col_accumulator[i] += a_val * b_val;
      }
    }

    std::ranges::sort(non_zero_indices);

    for (int i : non_zero_indices) {
      if (std::abs(col_accumulator[i]) > 1e-9) {
        thread_val.push_back(col_accumulator[i]);
        thread_row_idx.push_back(i);
        current_nnz++;
      }
      col_accumulator[i] = {0.0, 0.0};
      is_non_zero[i] = false;
    }
    non_zero_indices.clear();

    thread_col_ptr[j - start_col + 1] = current_nnz;
  }
}

}  // namespace

BorunovVComplexCcsSTL::BorunovVComplexCcsSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().resize(1);
}

bool BorunovVComplexCcsSTL::ValidationImpl() {
  const auto &a = GetInput().first;
  const auto &b = GetInput().second;
  if (a.num_cols != b.num_rows) {
    return false;
  }
  if (a.col_ptrs.size() != static_cast<size_t>(a.num_cols) + 1 ||
      b.col_ptrs.size() != static_cast<size_t>(b.num_cols) + 1) {
    return false;
  }
  return true;
}

bool BorunovVComplexCcsSTL::PreProcessingImpl() {
  const auto &a = GetInput().first;
  const auto &b = GetInput().second;
  auto &c = GetOutput()[0];

  c.num_rows = a.num_rows;
  c.num_cols = b.num_cols;
  c.col_ptrs.assign(c.num_cols + 1, 0);
  c.values.clear();
  c.row_indices.clear();

  return true;
}

bool BorunovVComplexCcsSTL::RunImpl() {
  const auto &a = GetInput().first;
  const auto &b = GetInput().second;
  auto &c = GetOutput()[0];

  int num_cols = b.num_cols;
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 4;
  }

  if (std::cmp_greater(num_threads, num_cols)) {
    num_threads = num_cols;
  }

  if (num_cols == 0) {
    return true;
  }

  std::vector<std::vector<std::complex<double>>> thread_values(num_threads);
  std::vector<std::vector<int>> thread_row_indices(num_threads);
  std::vector<std::vector<int>> thread_col_ptrs(num_threads);

  auto worker = [&](int thread_id) {
    WorkerThread(thread_id, static_cast<int>(num_threads), num_cols, a, b, thread_values[thread_id],
                 thread_row_indices[thread_id], thread_col_ptrs[thread_id]);
  };

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  for (unsigned int i = 0; i < num_threads; ++i) {
    threads.emplace_back(worker, i);
  }

  for (auto &t : threads) {
    t.join();
  }

  // merge results
  int total_nnz = 0;
  for (unsigned int i = 0; i < num_threads; ++i) {
    total_nnz += static_cast<int>(thread_values[i].size());
  }

  c.values.reserve(total_nnz);
  c.row_indices.reserve(total_nnz);

  int current_global_ptr = 0;
  int num_threads_int = static_cast<int>(num_threads);
  for (unsigned int i = 0; i < num_threads; ++i) {
    int i_int = static_cast<int>(i);
    int start_col = (num_cols * i_int) / num_threads_int;
    int end_col = (num_cols * (i_int + 1)) / num_threads_int;

    for (int j = 0; j < end_col - start_col; ++j) {
      c.col_ptrs[start_col + j + 1] = current_global_ptr + thread_col_ptrs[i][j + 1];
    }

    current_global_ptr += static_cast<int>(thread_values[i].size());

    c.values.insert(c.values.end(), thread_values[i].begin(), thread_values[i].end());
    c.row_indices.insert(c.row_indices.end(), thread_row_indices[i].begin(), thread_row_indices[i].end());
  }

  return true;
}

bool BorunovVComplexCcsSTL::PostProcessingImpl() {
  return true;
}

}  // namespace borunov_v_complex_ccs
