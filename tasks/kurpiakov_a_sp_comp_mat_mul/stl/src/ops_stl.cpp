#include "kurpiakov_a_sp_comp_mat_mul/stl/include/ops_stl.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include "kurpiakov_a_sp_comp_mat_mul/common/include/common.hpp"
#include "util/include/util.hpp"

namespace kurpiakov_a_sp_comp_mat_mul {

namespace {

bool ValidateCSR(const SparseMatrix &m) {
  if (m.rows <= 0 || m.cols <= 0) {
    return false;
  }
  if (static_cast<int>(m.row_ptr.size()) != m.rows + 1) {
    return false;
  }
  if (m.row_ptr[0] != 0) {
    return false;
  }
  if (std::cmp_not_equal(m.values.size(), m.row_ptr[m.rows])) {
    return false;
  }
  if (m.col_indices.size() != m.values.size()) {
    return false;
  }
  for (int i = 0; i < m.rows; ++i) {
    for (int j = m.row_ptr[i]; j < m.row_ptr[i + 1]; ++j) {
      if (m.col_indices[j] < 0 || m.col_indices[j] >= m.cols) {
        return false;
      }
    }
  }
  return true;
}

struct ThreadLocalRowState {
  explicit ThreadLocalRowState(int cols) : row_acc(cols), row_used(cols, 0) {}

  std::vector<ComplexD> row_acc;
  std::vector<char> row_used;
  std::vector<int> used_cols;
};

void MultiplySingleRow(const SparseMatrix &a, const SparseMatrix &b, int row_idx, ThreadLocalRowState &state,
                       std::vector<ComplexD> &out_values, std::vector<int> &out_cols) {
  state.used_cols.clear();

  for (int ja = a.row_ptr[row_idx]; ja < a.row_ptr[row_idx + 1]; ++ja) {
    const int ka = a.col_indices[ja];
    const ComplexD &a_val = a.values[ja];

    for (int jb = b.row_ptr[ka]; jb < b.row_ptr[ka + 1]; ++jb) {
      const int cb = b.col_indices[jb];
      const ComplexD &b_val = b.values[jb];

      if (state.row_used[cb] == 0) {
        state.row_used[cb] = 1;
        state.row_acc[cb] = ComplexD();
        state.used_cols.push_back(cb);
      }

      state.row_acc[cb] += a_val * b_val;
    }
  }

  std::ranges::sort(state.used_cols);

  out_values.clear();
  out_cols.clear();
  out_values.reserve(state.used_cols.size());
  out_cols.reserve(state.used_cols.size());

  for (int col : state.used_cols) {
    out_values.push_back(state.row_acc[col]);
    out_cols.push_back(col);
    state.row_used[col] = 0;
  }
}

void ComputeRowsWithThreads(const SparseMatrix &a, const SparseMatrix &b, int rows, int num_threads,
                            std::vector<std::vector<ComplexD>> &row_values, std::vector<std::vector<int>> &row_cols) {
  std::atomic<int> next_row(0);
  std::vector<std::thread> workers;
  workers.reserve(num_threads);

  for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    workers.emplace_back([&]() {
      ThreadLocalRowState state(b.cols);

      while (true) {
        const int row_idx = next_row.fetch_add(1, std::memory_order_relaxed);
        if (row_idx >= rows) {
          break;
        }

        MultiplySingleRow(a, b, row_idx, state, row_values[row_idx], row_cols[row_idx]);
      }
    });
  }

  for (auto &worker : workers) {
    worker.join();
  }
}

SparseMatrix BuildResult(int rows, int cols, const std::vector<std::vector<ComplexD>> &row_values,
                         const std::vector<std::vector<int>> &row_cols) {
  SparseMatrix result(rows, cols);

  std::size_t total_nnz = 0;
  for (int i = 0; i < rows; ++i) {
    total_nnz += row_values[i].size();
  }

  result.values.reserve(total_nnz);
  result.col_indices.reserve(total_nnz);

  for (int i = 0; i < rows; ++i) {
    result.values.insert(result.values.end(), row_values[i].begin(), row_values[i].end());
    result.col_indices.insert(result.col_indices.end(), row_cols[i].begin(), row_cols[i].end());
    result.row_ptr[i + 1] = static_cast<int>(result.values.size());
  }

  return result;
}

}  // namespace

KurpiakovACRSMatMulSTL::KurpiakovACRSMatMulSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = SparseMatrix();
}

bool KurpiakovACRSMatMulSTL::ValidationImpl() {
  const auto &[a, b] = GetInput();

  if (!ValidateCSR(a) || !ValidateCSR(b)) {
    return false;
  }

  return a.cols == b.rows;
}

bool KurpiakovACRSMatMulSTL::PreProcessingImpl() {
  return true;
}

bool KurpiakovACRSMatMulSTL::RunImpl() {
  const auto &[a, b] = GetInput();
  const int rows = a.rows;
  const int cols = b.cols;

  const int requested_threads = ppc::util::GetNumThreads();
  const int num_threads = std::max(1, std::min(requested_threads, rows));

  std::vector<std::vector<ComplexD>> row_values(rows);
  std::vector<std::vector<int>> row_cols(rows);

  ComputeRowsWithThreads(a, b, rows, num_threads, row_values, row_cols);
  GetOutput() = BuildResult(rows, cols, row_values, row_cols);
  return true;
}

bool KurpiakovACRSMatMulSTL::PostProcessingImpl() {
  return true;
}

}  // namespace kurpiakov_a_sp_comp_mat_mul
