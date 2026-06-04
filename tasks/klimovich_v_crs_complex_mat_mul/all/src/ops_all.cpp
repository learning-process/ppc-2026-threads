#include "klimovich_v_crs_complex_mat_mul/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "klimovich_v_crs_complex_mat_mul/common/include/common.hpp"

namespace klimovich_v_crs_complex_mat_mul {
namespace {

struct RowStage {
  std::vector<int> cols;
  std::vector<Cplx> vals;
};

void GustavsonRow(const CrsMatrix &lhs, const CrsMatrix &rhs, int row, std::vector<Cplx> &spa,
                  std::vector<int> &touched_by_row, std::vector<int> &touched_cols, RowStage &stage) {
  touched_cols.clear();

  for (int lp = lhs.row_offsets[row]; lp < lhs.row_offsets[row + 1]; ++lp) {
    const int k = lhs.col_indices[lp];
    const Cplx a_ik = lhs.data[lp];
    for (int rq = rhs.row_offsets[k]; rq < rhs.row_offsets[k + 1]; ++rq) {
      const int j = rhs.col_indices[rq];
      if (touched_by_row[j] != row) {
        touched_by_row[j] = row;
        touched_cols.push_back(j);
        spa[j] = a_ik * rhs.data[rq];
      } else {
        spa[j] += a_ik * rhs.data[rq];
      }
    }
  }

  std::ranges::sort(touched_cols);

  stage.cols.clear();
  stage.vals.clear();
  stage.cols.reserve(touched_cols.size());
  stage.vals.reserve(touched_cols.size());

  for (const int j : touched_cols) {
    const Cplx v = spa[j];
    spa[j] = Cplx(0.0, 0.0);
    if (std::abs(v.real()) > kZeroDropTol || std::abs(v.imag()) > kZeroDropTol) {
      stage.cols.push_back(j);
      stage.vals.push_back(v);
    }
  }
}

void RowRange(int total_rows, int world_size, int rank, int &begin, int &end) {
  const int base = total_rows / world_size;
  const int extra = total_rows % world_size;
  begin = (rank * base) + std::min(rank, extra);
  end = begin + base + (rank < extra ? 1 : 0);
}

// Helpers to reduce cognitive complexity of RunImpl
void FillRowsPerProc(int lhs_n_rows, int world_size, int rank, std::vector<int> &rows_per_proc,
                     std::vector<int> &rows_displs) {
  if (rank == 0) {
    rows_per_proc.assign(static_cast<std::size_t>(world_size), 0);
    rows_displs.assign(static_cast<std::size_t>(world_size), 0);
    for (int proc = 0; proc < world_size; ++proc) {
      int b = 0;
      int e = 0;
      RowRange(lhs_n_rows, world_size, proc, b, e);
      rows_per_proc[proc] = e - b;
      rows_displs[proc] = b;
    }
  }
}

int GatherPayloadCountsAndDispls(int local_payload, int world_size, int rank, std::vector<int> &payload_counts,
                                 std::vector<int> &payload_displs) {
  if (rank == 0) {
    payload_counts.assign(static_cast<std::size_t>(world_size), 0);
    payload_displs.assign(static_cast<std::size_t>(world_size), 0);
  }
  MPI_Gather(&local_payload, 1, MPI_INT, rank == 0 ? payload_counts.data() : nullptr, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int total_payload = 0;
  if (rank == 0) {
    for (int proc = 0; proc < world_size; ++proc) {
      payload_displs[proc] = total_payload;
      total_payload += payload_counts[proc];
    }
  }
  return total_payload;
}

void BuildPayloadCountsD(const std::vector<int> &payload_counts, const std::vector<int> &payload_displs, int world_size,
                         int rank, std::vector<int> &payload_counts_d, std::vector<int> &payload_displs_d) {
  if (rank == 0) {
    payload_counts_d.assign(static_cast<std::size_t>(world_size), 0);
    payload_displs_d.assign(static_cast<std::size_t>(world_size), 0);
    for (int proc = 0; proc < world_size; ++proc) {
      payload_counts_d[proc] = payload_counts[proc] * 2;
      payload_displs_d[proc] = payload_displs[proc] * 2;
    }
  }
}

}  // namespace

void KlimovichVCrsComplexMatMulAll::BroadcastOperand(CrsMatrix &m, int root) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::array<int, 3> meta{0, 0, 0};
  if (rank == root) {
    meta[0] = m.n_rows;
    meta[1] = m.n_cols;
    meta[2] = static_cast<int>(m.data.size());
  }
  MPI_Bcast(meta.data(), 3, MPI_INT, root, MPI_COMM_WORLD);

  if (rank != root) {
    m.n_rows = meta[0];
    m.n_cols = meta[1];
    m.row_offsets.assign(static_cast<std::size_t>(meta[0]) + 1, 0);
    m.col_indices.assign(static_cast<std::size_t>(meta[2]), 0);
    m.data.assign(static_cast<std::size_t>(meta[2]), Cplx(0.0, 0.0));
  }

  if (meta[0] > 0) {
    MPI_Bcast(m.row_offsets.data(), meta[0] + 1, MPI_INT, root, MPI_COMM_WORLD);
  }
  if (meta[2] > 0) {
    MPI_Bcast(m.col_indices.data(), meta[2], MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(m.data.data(), meta[2] * 2, MPI_DOUBLE, root, MPI_COMM_WORLD);
  }
}

void KlimovichVCrsComplexMatMulAll::ComputeLocalRows(const CrsMatrix &lhs, const CrsMatrix &rhs, int row_begin,
                                                     int row_end, std::vector<int> &local_nnz_per_row,
                                                     std::vector<int> &local_cols, std::vector<Cplx> &local_vals) {
  const int local_rows = row_end - row_begin;
  std::vector<RowStage> stages(static_cast<std::size_t>(local_rows));

#pragma omp parallel default(none) shared(lhs, rhs, stages, row_begin, row_end)
  {
    std::vector<Cplx> spa(static_cast<std::size_t>(rhs.n_cols));
    std::vector<int> touched_by_row(static_cast<std::size_t>(rhs.n_cols), -1);
    std::vector<int> touched_cols;
    touched_cols.reserve(static_cast<std::size_t>(rhs.n_cols));

#pragma omp for schedule(dynamic, 16)
    for (int i = row_begin; i < row_end; ++i) {
      GustavsonRow(lhs, rhs, i, spa, touched_by_row, touched_cols, stages[i - row_begin]);
    }
  }

  local_nnz_per_row.assign(static_cast<std::size_t>(local_rows), 0);
  std::size_t total = 0;
  for (int i = 0; i < local_rows; ++i) {
    local_nnz_per_row[i] = static_cast<int>(stages[i].cols.size());
    total += stages[i].cols.size();
  }

  local_cols.clear();
  local_vals.clear();
  local_cols.reserve(total);
  local_vals.reserve(total);
  for (auto &stage : stages) {
    local_cols.insert(local_cols.end(), stage.cols.begin(), stage.cols.end());
    local_vals.insert(local_vals.end(), stage.vals.begin(), stage.vals.end());
  }
}

KlimovichVCrsComplexMatMulAll::KlimovichVCrsComplexMatMulAll(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    GetInput() = in;
  }
  GetOutput() = CrsMatrix();
}

bool KlimovichVCrsComplexMatMulAll::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) {
    return true;
  }
  const auto &lhs = std::get<0>(GetInput());
  const auto &rhs = std::get<1>(GetInput());
  return lhs.n_cols == rhs.n_rows;
}

bool KlimovichVCrsComplexMatMulAll::PreProcessingImpl() {
  return true;
}

bool KlimovichVCrsComplexMatMulAll::RunImpl() {
  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  CrsMatrix lhs;
  CrsMatrix rhs;
  if (rank == 0) {
    lhs = std::get<0>(GetInput());
    rhs = std::get<1>(GetInput());
  }
  BroadcastOperand(lhs, 0);
  BroadcastOperand(rhs, 0);

  int row_begin = 0;
  int row_end = 0;
  RowRange(lhs.n_rows, world_size, rank, row_begin, row_end);

  std::vector<int> local_nnz_per_row;
  std::vector<int> local_cols;
  std::vector<Cplx> local_vals;
  ComputeLocalRows(lhs, rhs, row_begin, row_end, local_nnz_per_row, local_cols, local_vals);

  std::vector<int> rows_per_proc;
  std::vector<int> rows_displs;
  FillRowsPerProc(lhs.n_rows, world_size, rank, rows_per_proc, rows_displs);

  std::vector<int> global_nnz_per_row;
  if (rank == 0) {
    global_nnz_per_row.assign(static_cast<std::size_t>(lhs.n_rows), 0);
  }
  MPI_Gatherv(local_nnz_per_row.data(), static_cast<int>(local_nnz_per_row.size()), MPI_INT,
              rank == 0 ? global_nnz_per_row.data() : nullptr, rank == 0 ? rows_per_proc.data() : nullptr,
              rank == 0 ? rows_displs.data() : nullptr, MPI_INT, 0, MPI_COMM_WORLD);

  const int local_payload = static_cast<int>(local_cols.size());
  std::vector<int> payload_counts;
  std::vector<int> payload_displs;
  int total_payload = GatherPayloadCountsAndDispls(local_payload, world_size, rank, payload_counts, payload_displs);

  std::vector<int> all_cols;
  std::vector<Cplx> all_vals;
  if (rank == 0) {
    all_cols.assign(static_cast<std::size_t>(total_payload), 0);
    all_vals.assign(static_cast<std::size_t>(total_payload), Cplx(0.0, 0.0));
  }

  MPI_Gatherv(local_cols.data(), local_payload, MPI_INT, rank == 0 ? all_cols.data() : nullptr,
              rank == 0 ? payload_counts.data() : nullptr, rank == 0 ? payload_displs.data() : nullptr, MPI_INT, 0,
              MPI_COMM_WORLD);

  std::vector<int> payload_counts_d;
  std::vector<int> payload_displs_d;
  BuildPayloadCountsD(payload_counts, payload_displs, world_size, rank, payload_counts_d, payload_displs_d);
  MPI_Gatherv(local_vals.data(), local_payload * 2, MPI_DOUBLE, rank == 0 ? all_vals.data() : nullptr,
              rank == 0 ? payload_counts_d.data() : nullptr, rank == 0 ? payload_displs_d.data() : nullptr, MPI_DOUBLE,
              0, MPI_COMM_WORLD);

  if (rank == 0) {
    CrsMatrix &out = GetOutput();
    out = CrsMatrix(lhs.n_rows, rhs.n_cols);
    for (int i = 0; i < lhs.n_rows; ++i) {
      out.row_offsets[i + 1] = out.row_offsets[i] + global_nnz_per_row[i];
    }
    out.col_indices = std::move(all_cols);
    out.data = std::move(all_vals);
  }

  return true;
}

bool KlimovichVCrsComplexMatMulAll::PostProcessingImpl() {
  return true;
}

}  // namespace klimovich_v_crs_complex_mat_mul
