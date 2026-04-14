#include "kurpiakov_a_sp_comp_mat_mul/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <atomic>
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

std::pair<int, int> GetRowRange(int total_rows, int rank, int size) {
  const int begin = (total_rows * rank) / size;
  const int end = (total_rows * (rank + 1)) / size;
  return {begin, end};
}

void MultiplySingleRow(const SparseMatrix &a, const SparseMatrix &b, int row_idx, std::vector<ComplexD> &row_acc,
                       std::vector<char> &row_used, std::vector<int> &used_cols, std::vector<ComplexD> &out_values,
                       std::vector<int> &out_cols) {
  used_cols.clear();

  for (int ja = a.row_ptr[row_idx]; ja < a.row_ptr[row_idx + 1]; ++ja) {
    const int ka = a.col_indices[ja];
    const ComplexD &a_val = a.values[ja];

    for (int jb = b.row_ptr[ka]; jb < b.row_ptr[ka + 1]; ++jb) {
      const int cb = b.col_indices[jb];
      const ComplexD &b_val = b.values[jb];

      if (row_used[cb] == 0) {
        row_used[cb] = 1;
        row_acc[cb] = ComplexD();
        used_cols.push_back(cb);
      }

      row_acc[cb] += a_val * b_val;
    }
  }

  std::sort(used_cols.begin(), used_cols.end());

  out_values.clear();
  out_cols.clear();
  out_values.reserve(used_cols.size());
  out_cols.reserve(used_cols.size());

  for (int col : used_cols) {
    out_values.push_back(row_acc[col]);
    out_cols.push_back(col);
    row_used[col] = 0;
  }
}

void ComputeLocalRowsThreads(const SparseMatrix &a, const SparseMatrix &b, int row_begin, int row_end,
                             std::vector<std::vector<ComplexD>> &local_values,
                             std::vector<std::vector<int>> &local_cols) {
  const int local_rows = row_end - row_begin;
  const int requested_threads = ppc::util::GetNumThreads();
  const int max_threads = std::max(1, local_rows);
  const int num_threads = std::max(1, std::min(requested_threads, max_threads));

  std::atomic<int> next_row(row_begin);
  std::vector<std::thread> workers;
  workers.reserve(num_threads);

  for (int t = 0; t < num_threads; ++t) {
    workers.emplace_back([&]() {
      std::vector<ComplexD> row_acc(b.cols);
      std::vector<char> row_used(b.cols, 0);
      std::vector<int> used_cols;

      while (true) {
        const int row = next_row.fetch_add(1, std::memory_order_relaxed);
        if (row >= row_end) {
          break;
        }

        const int local_idx = row - row_begin;
        MultiplySingleRow(a, b, row, row_acc, row_used, used_cols, local_values[local_idx], local_cols[local_idx]);
      }
    });
  }

  for (auto &worker : workers) {
    worker.join();
  }
}

}  // namespace

KurpiakovACRSMatMulALL::KurpiakovACRSMatMulALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = SparseMatrix();
}

bool KurpiakovACRSMatMulALL::ValidationImpl() {
  const auto &[a, b] = GetInput();

  if (!ValidateCSR(a) || !ValidateCSR(b)) {
    return false;
  }

  return a.cols == b.rows;
}

bool KurpiakovACRSMatMulALL::PreProcessingImpl() {
  return true;
}

bool KurpiakovACRSMatMulALL::RunImpl() {
  const auto &[a, b] = GetInput();
  const int rows = a.rows;
  const int cols = b.cols;

  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  const auto [row_begin, row_end] = GetRowRange(rows, rank, world_size);
  const int local_rows = row_end - row_begin;

  std::vector<std::vector<ComplexD>> local_values(local_rows);
  std::vector<std::vector<int>> local_cols(local_rows);

  ComputeLocalRowsThreads(a, b, row_begin, row_end, local_values, local_cols);

  std::vector<int> local_row_nnz(rows, 0);
  for (int local_i = 0; local_i < local_rows; ++local_i) {
    local_row_nnz[row_begin + local_i] = static_cast<int>(local_values[local_i].size());
  }

  std::vector<int> global_row_nnz(rows, 0);
  MPI_Allreduce(local_row_nnz.data(), global_row_nnz.data(), rows, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  std::vector<int> global_row_ptr(rows + 1, 0);
  for (int i = 0; i < rows; ++i) {
    global_row_ptr[i + 1] = global_row_ptr[i] + global_row_nnz[i];
  }

  const int total_nnz = global_row_ptr[rows];
  const int local_nnz = global_row_ptr[row_end] - global_row_ptr[row_begin];

  std::vector<double> local_re(local_nnz);
  std::vector<double> local_im(local_nnz);
  std::vector<int> local_col_indices(local_nnz);

  int pos = 0;
  for (int local_i = 0; local_i < local_rows; ++local_i) {
    const auto &vals = local_values[local_i];
    const auto &cols_row = local_cols[local_i];

    for (std::size_t j = 0; j < vals.size(); ++j) {
      local_re[pos] = vals[j].re;
      local_im[pos] = vals[j].im;
      local_col_indices[pos] = cols_row[j];
      ++pos;
    }
  }

  std::vector<int> recv_counts(world_size, 0);
  MPI_Allgather(&local_nnz, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> recv_displs(world_size, 0);
  for (int r = 1; r < world_size; ++r) {
    recv_displs[r] = recv_displs[r - 1] + recv_counts[r - 1];
  }

  std::vector<double> global_re;
  std::vector<double> global_im;
  std::vector<int> global_col_indices;

  if (rank == 0) {
    global_re.resize(total_nnz);
    global_im.resize(total_nnz);
    global_col_indices.resize(total_nnz);
  }

  MPI_Gatherv(local_re.data(), local_nnz, MPI_DOUBLE, global_re.data(), recv_counts.data(), recv_displs.data(),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(local_im.data(), local_nnz, MPI_DOUBLE, global_im.data(), recv_counts.data(), recv_displs.data(),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(local_col_indices.data(), local_nnz, MPI_INT, global_col_indices.data(), recv_counts.data(),
              recv_displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    global_re.resize(total_nnz);
    global_im.resize(total_nnz);
    global_col_indices.resize(total_nnz);
  }

  MPI_Bcast(global_re.data(), total_nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(global_im.data(), total_nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(global_col_indices.data(), total_nnz, MPI_INT, 0, MPI_COMM_WORLD);

  SparseMatrix result(rows, cols);
  result.row_ptr = std::move(global_row_ptr);
  result.col_indices = std::move(global_col_indices);
  result.values.resize(static_cast<std::size_t>(total_nnz));

  for (int i = 0; i < total_nnz; ++i) {
    result.values[static_cast<std::size_t>(i)] = ComplexD(global_re[i], global_im[i]);
  }

  GetOutput() = std::move(result);
  return true;
}

bool KurpiakovACRSMatMulALL::PostProcessingImpl() {
  return true;
}

}  // namespace kurpiakov_a_sp_comp_mat_mul
