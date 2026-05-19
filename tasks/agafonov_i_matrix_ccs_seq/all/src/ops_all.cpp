#include "agafonov_i_matrix_ccs_seq/all/include/ops_all.hpp"

#include <mpi.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "agafonov_i_matrix_ccs_seq/common/include/common.hpp"

namespace agafonov_i_matrix_ccs_seq {

namespace {

void AssembleFinalMatrix(InType::first_type &c, size_t total_cols, size_t chunk, size_t rem, int size,
                         const std::vector<int> &all_col_nnz, const std::vector<int> &displs,
                         const std::vector<double> &all_vals, const std::vector<int> &all_rind) {
  std::vector<size_t> proc_col_start(size);
  std::vector<size_t> proc_col_end(size);
  for (int proc_idx = 0; proc_idx < size; ++proc_idx) {
    proc_col_start[proc_idx] = (static_cast<size_t>(proc_idx) * chunk) + std::min(static_cast<size_t>(proc_idx), rem);
    proc_col_end[proc_idx] = proc_col_start[proc_idx] + chunk + (std::cmp_less(proc_idx, rem) ? 1 : 0);
  }

  std::vector<std::vector<double>> col_vals(total_cols);
  std::vector<std::vector<int>> col_rind(total_cols);
  for (int proc_idx = 0; proc_idx < size; ++proc_idx) {
    int offset = displs[proc_idx];
    for (size_t j = proc_col_start[proc_idx]; j < proc_col_end[proc_idx]; ++j) {
      int n = all_col_nnz[j];
      col_vals[j].assign(all_vals.begin() + offset, all_vals.begin() + offset + n);
      col_rind[j].assign(all_rind.begin() + offset, all_rind.begin() + offset + n);
      offset += n;
    }
  }

  int cur_nnz = 0;
  for (size_t j = 0; j < total_cols; ++j) {
    c.col_ptrs[j] = cur_nnz;
    c.vals.insert(c.vals.end(), col_vals[j].begin(), col_vals[j].end());
    c.row_inds.insert(c.row_inds.end(), col_rind[j].begin(), col_rind[j].end());
    cur_nnz += static_cast<int>(col_vals[j].size());
  }
  c.col_ptrs[total_cols] = cur_nnz;
  c.nnz = cur_nnz;
}

}  // namespace

AgafonovIMatrixCCSALL::AgafonovIMatrixCCSALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool AgafonovIMatrixCCSALL::ValidationImpl() {
  const auto &left = GetInput().first;
  const auto &right = GetInput().second;
  return (left.cols_num == right.rows_num) && (left.col_ptrs.size() == left.cols_num + 1) &&
         (right.col_ptrs.size() == right.cols_num + 1);
}

bool AgafonovIMatrixCCSALL::PreProcessingImpl() {
  GetOutput().vals.clear();
  GetOutput().row_inds.clear();
  return true;
}

void AgafonovIMatrixCCSALL::ProcessColumn(size_t j, const InType::first_type &a, const InType::second_type &b,
                                          std::vector<double> &accumulator, std::vector<size_t> &active_rows,
                                          std::vector<bool> &row_mask, std::vector<double> &local_v,
                                          std::vector<int> &local_r) {
  const auto b_col_start = static_cast<size_t>(b.col_ptrs[j]);
  const auto b_col_end = static_cast<size_t>(b.col_ptrs[j + 1]);
  if (b_col_start == b_col_end) {
    return;
  }

  for (size_t kb = b_col_start; kb < b_col_end; ++kb) {
    const auto k = static_cast<size_t>(b.row_inds[kb]);
    const double v_b = b.vals[kb];
    const auto a_col_start = static_cast<size_t>(a.col_ptrs[k]);
    const auto a_col_end = static_cast<size_t>(a.col_ptrs[k + 1]);

    for (size_t ka = a_col_start; ka < a_col_end; ++ka) {
      const auto i = static_cast<size_t>(a.row_inds[ka]);
      if (!row_mask[i]) {
        row_mask[i] = true;
        active_rows.push_back(i);
      }
      accumulator[i] += a.vals[ka] * v_b;
    }
  }

  std::ranges::sort(active_rows);

  for (const auto row_idx : active_rows) {
    if (std::abs(accumulator[row_idx]) > 1e-15) {
      local_v.push_back(accumulator[row_idx]);
      local_r.push_back(static_cast<int>(row_idx));
    }
    row_mask[row_idx] = false;
    accumulator[row_idx] = 0.0;
  }
  active_rows.clear();
}

bool AgafonovIMatrixCCSALL::RunImpl() {
  const auto &a = GetInput().first;
  const auto &b = GetInput().second;
  auto &c = GetOutput();

  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  c.rows_num = a.rows_num;
  c.cols_num = b.cols_num;
  c.col_ptrs.assign(c.cols_num + 1, 0);

  const size_t total_cols = b.cols_num;
  const size_t chunk = total_cols / static_cast<size_t>(size);
  const size_t rem = total_cols % static_cast<size_t>(size);
  const size_t col_start = (static_cast<size_t>(rank) * chunk) + std::min(static_cast<size_t>(rank), rem);
  const size_t col_end = col_start + chunk + (std::cmp_less(rank, rem) ? 1 : 0);

  std::vector<std::vector<double>> local_vals(total_cols);
  std::vector<std::vector<int>> local_rows(total_cols);

  tbb::parallel_for(tbb::blocked_range<size_t>(col_start, col_end), [&](const tbb::blocked_range<size_t> &range) {
    std::vector<double> accumulator(a.rows_num, 0.0);
    std::vector<size_t> active_rows;
    std::vector<bool> row_mask(a.rows_num, false);

    for (size_t j = range.begin(); j < range.end(); ++j) {
      ProcessColumn(j, a, b, accumulator, active_rows, row_mask, local_vals[j], local_rows[j]);
    }
  });

  std::vector<int> my_col_nnz(total_cols, 0);
  for (size_t j = col_start; j < col_end; ++j) {
    my_col_nnz[j] = static_cast<int>(local_vals[j].size());
  }

  std::vector<int> all_col_nnz;
  if (rank == 0) {
    all_col_nnz.resize(total_cols, 0);
  }

  MPI_Reduce(my_col_nnz.data(), all_col_nnz.data(), static_cast<int>(total_cols), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  std::vector<double> my_vals;
  std::vector<int> my_rind;
  for (size_t j = col_start; j < col_end; ++j) {
    my_vals.insert(my_vals.end(), local_vals[j].begin(), local_vals[j].end());
    my_rind.insert(my_rind.end(), local_rows[j].begin(), local_rows[j].end());
  }

  int my_count = static_cast<int>(my_vals.size());
  std::vector<int> counts(size);
  MPI_Gather(&my_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> displs(size, 0);
  int total_nnz = 0;
  if (rank == 0) {
    for (int i = 0; i < size; ++i) {
      displs[i] = total_nnz;
      total_nnz += counts[i];
    }
  }

  std::vector<double> all_vals;
  std::vector<int> all_rind;
  if (rank == 0) {
    all_vals.resize(total_nnz);
    all_rind.resize(total_nnz);
  }

  MPI_Gatherv(my_vals.data(), my_count, MPI_DOUBLE, all_vals.data(), counts.data(), displs.data(), MPI_DOUBLE, 0,
              MPI_COMM_WORLD);
  MPI_Gatherv(my_rind.data(), my_count, MPI_INT, all_rind.data(), counts.data(), displs.data(), MPI_INT, 0,
              MPI_COMM_WORLD);

  if (rank == 0) {
    AssembleFinalMatrix(c, total_cols, chunk, rem, size, all_col_nnz, displs, all_vals, all_rind);
  }

  MPI_Bcast(&c.rows_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&c.cols_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&c.nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    c.col_ptrs.resize(c.cols_num + 1);
    c.vals.resize(c.nnz);
    c.row_inds.resize(c.nnz);
  }

  MPI_Bcast(c.col_ptrs.data(), static_cast<int>(c.col_ptrs.size()), MPI_INT, 0, MPI_COMM_WORLD);

  if (c.nnz > 0) {
    MPI_Bcast(c.vals.data(), static_cast<int>(c.nnz), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(c.row_inds.data(), static_cast<int>(c.nnz), MPI_INT, 0, MPI_COMM_WORLD);
  }
  return true;
}

bool AgafonovIMatrixCCSALL::PostProcessingImpl() {
  return true;
}

}  // namespace agafonov_i_matrix_ccs_seq
