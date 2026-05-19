#include "kapanova_s_sparse_matrix_mult_ccs/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "kapanova_s_sparse_matrix_mult_ccs/common/include/common.hpp"

namespace kapanova_s_sparse_matrix_mult_ccs {

KapanovaSSparseMatrixMultCCSALL::KapanovaSSparseMatrixMultCCSALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KapanovaSSparseMatrixMultCCSALL::ValidationImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  return (a.cols == b.rows && a.rows > 0 && a.cols > 0 && b.rows > 0 && b.cols > 0 &&
          a.col_ptrs.size() == static_cast<size_t>(a.cols + 1) && b.col_ptrs.size() == static_cast<size_t>(b.cols + 1));
}

bool KapanovaSSparseMatrixMultCCSALL::PreProcessingImpl() {
  return true;
}
bool KapanovaSSparseMatrixMultCCSALL::PostProcessingImpl() {
  return true;
}

namespace {

using MpiU64 = std::uint64_t;
MPI_Datatype k_mpi_u64 = MPI_UINT64_T;

std::vector<MpiU64> ComputeBalancedRanges(int total_cols, int num_procs, const CCSMatrix &a, const CCSMatrix &b) {
  std::vector<MpiU64> ranges(static_cast<size_t>(num_procs) + 1, 0);
  ranges[num_procs] = static_cast<MpiU64>(total_cols);
  if (total_cols == 0) {
    return ranges;
  }

  std::vector<int> cost(static_cast<size_t>(total_cols), 0);
  size_t total_cost = 0;
#pragma omp parallel for reduction(+ : total_cost) schedule(guided) default(none) shared(a, b, cost, total_cols)
  for (int col = 0; col < total_cols; ++col) {
    int c = 0;
    for (size_t k = b.col_ptrs[col]; k < b.col_ptrs[col + 1]; ++k) {
      c += static_cast<int>(a.col_ptrs[b.row_indices[k] + 1] - a.col_ptrs[b.row_indices[k]]);
    }
    cost[static_cast<size_t>(col)] = c;
    total_cost += static_cast<size_t>(c);
  }
  if (total_cost == 0) {
    return ranges;
  }
  size_t per = total_cost / static_cast<size_t>(num_procs);
  size_t cur = 0;
  size_t acc = 0;
  auto total_cols_sz = static_cast<size_t>(total_cols);
  for (int proc = 1; proc < num_procs; ++proc) {
    size_t target = static_cast<size_t>(proc) * per;
    while (cur < total_cols_sz && acc < target) {
      acc += static_cast<size_t>(cost[cur]);
      ++cur;
    }
    ranges[static_cast<size_t>(proc)] = static_cast<MpiU64>(cur);
  }
  return ranges;
}

void ProcessColumn(size_t gcol, const CCSMatrix &a, const CCSMatrix &b, std::vector<MpiU64> &out_rows,
                   std::vector<MpiU64> &out_cols, std::vector<double> &out_vals, double *accum, char *used,
                   size_t *active, int &active_count) {
  for (size_t k = b.col_ptrs[gcol]; k < b.col_ptrs[gcol + 1]; ++k) {
    size_t row_b = b.row_indices[k];
    double vb = b.values[k];
    for (size_t zc = a.col_ptrs[row_b]; zc < a.col_ptrs[row_b + 1]; ++zc) {
      size_t i = a.row_indices[zc];
      double va = a.values[zc];
      if (used[i] == 0) {
        used[i] = 1;
        active[active_count++] = i;
        accum[i] = va * vb;
      } else {
        accum[i] += va * vb;
      }
    }
  }
  for (int idx = 0; idx < active_count; ++idx) {
    size_t i = active[static_cast<size_t>(idx)];
    if (accum[i] != 0.0) {
      out_rows.push_back(static_cast<MpiU64>(i));
      out_cols.push_back(static_cast<MpiU64>(gcol));
      out_vals.push_back(accum[i]);
    }
    used[i] = 0;
    accum[i] = 0.0;
  }
}

void ComputeLocalColumns(size_t start, size_t local_cols, const CCSMatrix &a, const CCSMatrix &b,
                         std::vector<MpiU64> &send_rows, std::vector<MpiU64> &send_cols,
                         std::vector<double> &send_vals) {
#pragma omp parallel default(none) shared(a, b, start, local_cols, send_rows, send_cols, send_vals)
  {
    std::vector<double> accum(a.rows, 0.0);
    std::vector<char> used(static_cast<size_t>(a.rows), 0);
    std::vector<size_t> active(static_cast<size_t>(a.rows));
    std::vector<MpiU64> thr_rows;
    std::vector<MpiU64> thr_cols;
    std::vector<double> thr_vals;
    int active_count = 0;

#pragma omp for schedule(guided, 32) nowait
    for (size_t j = 0; j < local_cols; ++j) {
      active_count = 0;
      ProcessColumn(start + j, a, b, thr_rows, thr_cols, thr_vals, accum.data(), used.data(), active.data(),
                    active_count);
    }
#pragma omp critical
    {
      send_rows.insert(send_rows.end(), thr_rows.begin(), thr_rows.end());
      send_cols.insert(send_cols.end(), thr_cols.begin(), thr_cols.end());
      send_vals.insert(send_vals.end(), thr_vals.begin(), thr_vals.end());
    }
  }
}

void BuildCcsOnRoot(OutType &c, int total, std::vector<MpiU64> &recv_rows, std::vector<MpiU64> &recv_cols,
                    std::vector<double> &recv_vals) {
  c.nnz = static_cast<size_t>(total);
  c.col_ptrs.assign(static_cast<size_t>(c.cols) + 1, 0);
  c.row_indices.resize(static_cast<size_t>(total));
  c.values.resize(static_cast<size_t>(total));

  if (total > 0) {
    auto total_sz = static_cast<size_t>(total);
    std::vector<size_t> idx(total_sz);
    for (size_t i = 0; i < total_sz; ++i) {
      idx[i] = i;
    }
    std::ranges::sort(idx, [&](size_t a_idx, size_t b_idx) {
      auto lc = recv_cols[a_idx];
      auto rc = recv_cols[b_idx];
      if (lc != rc) {
        return lc < rc;
      }
      return recv_rows[a_idx] < recv_rows[b_idx];
    });
    for (size_t i = 0; i < total_sz; ++i) {
      size_t src = idx[i];
      c.row_indices[i] = static_cast<size_t>(recv_rows[src]);
      c.values[i] = recv_vals[src];
      c.col_ptrs[recv_cols[src] + 1]++;
    }
    for (size_t j = 0; j < static_cast<size_t>(c.cols); ++j) {
      c.col_ptrs[j + 1] += c.col_ptrs[j];
    }
  }
}

void BroadcastResult(OutType &c, int rank, MPI_Datatype mpi_type) {
  auto nnz_b = static_cast<MpiU64>(c.nnz);
  auto cols_b = static_cast<MpiU64>(c.cols);
  MPI_Bcast(&nnz_b, 1, mpi_type, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cols_b, 1, mpi_type, 0, MPI_COMM_WORLD);

  auto nnz_i = static_cast<int>(nnz_b);
  auto cols_i = static_cast<int>(cols_b);

  if (rank != 0) {
    c.nnz = static_cast<size_t>(nnz_b);
    c.cols = cols_i;
    c.col_ptrs.resize(static_cast<size_t>(cols_i) + 1);
    c.row_indices.resize(static_cast<size_t>(nnz_i));
    c.values.resize(static_cast<size_t>(nnz_i));
  }

  MPI_Bcast(c.col_ptrs.data(), cols_i + 1, mpi_type, 0, MPI_COMM_WORLD);
  MPI_Bcast(c.row_indices.data(), nnz_i, mpi_type, 0, MPI_COMM_WORLD);
  MPI_Bcast(c.values.data(), nnz_i, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

}  // namespace

bool KapanovaSSparseMatrixMultCCSALL::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  auto &c = GetOutput();
  c.rows = a.rows;
  c.cols = b.cols;

  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  auto r_sz = static_cast<size_t>(size);
  std::vector<MpiU64> ranges(r_sz + 1, 0);
  if (rank == 0) {
    ranges = ComputeBalancedRanges(static_cast<int>(c.cols), size, a, b);
  }
  MPI_Bcast(ranges.data(), size + 1, k_mpi_u64, 0, MPI_COMM_WORLD);

  auto start = static_cast<size_t>(ranges[rank]);
  auto local_cols = static_cast<size_t>(ranges[rank + 1]) - start;

  std::vector<MpiU64> send_rows;
  std::vector<MpiU64> send_cols;
  std::vector<double> send_vals;

  if (local_cols > 0) {
    ComputeLocalColumns(start, local_cols, a, b, send_rows, send_cols, send_vals);
  }

  int local_nnz = static_cast<int>(send_rows.size());
  std::vector<int> counts(r_sz, 0);
  MPI_Gather(&local_nnz, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> displs(r_sz, 0);
  int total = 0;
  if (rank == 0) {
    for (size_t i = 0; i < r_sz; ++i) {
      displs[i] = total;
      total += counts[i];
    }
  }

  auto total_sz = static_cast<size_t>(total);
  auto safe_sz = std::max(total_sz, static_cast<size_t>(1));
  std::vector<MpiU64> recv_rows(rank == 0 ? safe_sz : 0);
  std::vector<MpiU64> recv_cols(rank == 0 ? safe_sz : 0);
  std::vector<double> recv_vals(rank == 0 ? safe_sz : 0);

  const MpiU64 *rp = send_rows.empty() ? nullptr : send_rows.data();
  const MpiU64 *cp = send_cols.empty() ? nullptr : send_cols.data();
  const double *vp = send_vals.empty() ? nullptr : send_vals.data();

  MPI_Gatherv(rp, local_nnz, k_mpi_u64, recv_rows.data(), counts.data(), displs.data(), k_mpi_u64, 0, MPI_COMM_WORLD);
  MPI_Gatherv(cp, local_nnz, k_mpi_u64, recv_cols.data(), counts.data(), displs.data(), k_mpi_u64, 0, MPI_COMM_WORLD);
  MPI_Gatherv(vp, local_nnz, MPI_DOUBLE, recv_vals.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    BuildCcsOnRoot(c, total, recv_rows, recv_cols, recv_vals);
  }

  BroadcastResult(c, rank, k_mpi_u64);

  return true;
}

}  // namespace kapanova_s_sparse_matrix_mult_ccs
