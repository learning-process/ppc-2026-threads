#include "posternak_a_crs_mul_complex_matrix/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <unordered_map>
#include <utility>
#include <vector>

#include "posternak_a_crs_mul_complex_matrix/common/include/common.hpp"

namespace {

size_t ComputeRowNoZeroCount(const posternak_a_crs_mul_complex_matrix::CRSMatrix &a,
                             const posternak_a_crs_mul_complex_matrix::CRSMatrix &b, int row, double threshold) {
  std::unordered_map<int, std::complex<double>> row_sum;

  for (int idx_a = a.index_row[row]; idx_a < a.index_row[row + 1]; ++idx_a) {
    int col_a = a.index_col[idx_a];
    auto val_a = a.values[idx_a];

    for (int idx_b = b.index_row[col_a]; idx_b < b.index_row[col_a + 1]; ++idx_b) {
      int col_b = b.index_col[idx_b];
      auto val_b = b.values[idx_b];
      row_sum[col_b] += val_a * val_b;
    }
  }

  size_t local = 0;
  for (const auto &[col, val] : row_sum) {
    if (std::abs(val) > threshold) {
      ++local;
    }
  }
  return local;
}

void BuildResultStructure(posternak_a_crs_mul_complex_matrix::CRSMatrix &res, std::vector<size_t> &row_prefix) {
  for (int i = 1; i < res.rows; ++i) {
    row_prefix[i] += row_prefix[i - 1];
  }

  const size_t total = row_prefix.empty() ? 0 : row_prefix.back();
  res.values.resize(total);
  res.index_col.resize(total);
  res.index_row.resize(res.rows + 1);

  for (int i = 0; i <= res.rows; ++i) {
    res.index_row[i] = (i == 0 ? 0 : static_cast<int>(row_prefix[i - 1]));
  }
}

void ComputeAndWriteRow(const posternak_a_crs_mul_complex_matrix::CRSMatrix &a,
                        const posternak_a_crs_mul_complex_matrix::CRSMatrix &b,
                        posternak_a_crs_mul_complex_matrix::CRSMatrix &res, int row, double threshold) {
  std::unordered_map<int, std::complex<double>> row_sum;

  for (int idx_a = a.index_row[row]; idx_a < a.index_row[row + 1]; ++idx_a) {
    int col_a = a.index_col[idx_a];
    auto val_a = a.values[idx_a];

    for (int idx_b = b.index_row[col_a]; idx_b < b.index_row[col_a + 1]; ++idx_b) {
      int col_b = b.index_col[idx_b];
      auto val_b = b.values[idx_b];
      row_sum[col_b] += val_a * val_b;
    }
  }

  std::vector<std::pair<int, std::complex<double>>> sorted(row_sum.begin(), row_sum.end());

  std::ranges::sort(sorted, [](const auto &p1, const auto &p2) { return p1.first < p2.first; });

  size_t pos = res.index_row[row];
  for (const auto &[col_idx, value] : sorted) {
    if (std::abs(value) > threshold) {
      res.values[pos] = value;
      res.index_col[pos] = col_idx;
      ++pos;
    }
  }
}

void ComputeLocalCounts(const posternak_a_crs_mul_complex_matrix::CRSMatrix &a,
                        const posternak_a_crs_mul_complex_matrix::CRSMatrix &b, std::vector<size_t> &local_counts,
                        int local_start, int local_count, double threshold) {
#pragma omp parallel for schedule(dynamic) default(none) shared(a, b, local_counts, local_start, local_count, threshold)
  for (int i = 0; i < local_count; ++i) {
    local_counts[i] = ComputeRowNoZeroCount(a, b, local_start + i, threshold);
  }
}

void GatherCountsToRoot(const std::vector<size_t> &local_counts, int local_count, std::vector<size_t> &global_counts,
                        const std::vector<int> &recv_counts, const std::vector<int> &displs, int rank) {
  if (rank == 0) {
    std::vector<size_t> send_buf = local_counts;
    MPI_Gatherv(send_buf.data(), local_count, MPI_UNSIGNED_LONG_LONG, global_counts.data(), recv_counts.data(),
                displs.data(), MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
  } else {
    MPI_Gatherv(local_counts.data(), local_count, MPI_UNSIGNED_LONG_LONG, nullptr, nullptr, nullptr,
                MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
  }
}

void BuildAndBroadcastStructure(posternak_a_crs_mul_complex_matrix::CRSMatrix &res,
                                std::vector<size_t> &global_counts_copy, int rank) {
  if (rank == 0) {
    BuildResultStructure(res, global_counts_copy);
  }
  res.index_row.resize(res.rows + 1);
  MPI_Bcast(res.index_row.data(), static_cast<int>(res.index_row.size()), MPI_INT, 0, MPI_COMM_WORLD);
}

void ComputeLocalRows(const posternak_a_crs_mul_complex_matrix::CRSMatrix &a,
                      const posternak_a_crs_mul_complex_matrix::CRSMatrix &b,
                      posternak_a_crs_mul_complex_matrix::CRSMatrix &res, int local_start, int local_count,
                      double threshold) {
#pragma omp parallel for schedule(dynamic) default(none) shared(a, b, res, local_start, local_count, threshold)
  for (int i = 0; i < local_count; ++i) {
    ComputeAndWriteRow(a, b, res, local_start + i, threshold);
  }
}

void PrepareGatherParams(std::vector<int> &g_counts, std::vector<int> &g_displs,
                         const posternak_a_crs_mul_complex_matrix::CRSMatrix &res, int rows_per_proc, int rem,
                         int size) {
  for (int process = 0; process < size; ++process) {
    int process_start = (process * rows_per_proc) + std::min(process, rem);
    int process_end = process_start + rows_per_proc + (process < rem ? 1 : 0);
    g_displs[process] = res.index_row[process_start];
    g_counts[process] = res.index_row[process_end] - res.index_row[process_start];
  }
}

void GatherResultData(posternak_a_crs_mul_complex_matrix::CRSMatrix &res, int local_start, int local_nnz,
                      const std::vector<int> &g_counts, const std::vector<int> &g_displs, int rank) {
  if (rank == 0) {
    std::vector<std::complex<double>> local_values_copy(res.values.data() + res.index_row[local_start],
                                                        res.values.data() + res.index_row[local_start] + local_nnz);
    std::vector<int> local_index_copy(res.index_col.data() + res.index_row[local_start],
                                      res.index_col.data() + res.index_row[local_start] + local_nnz);

    MPI_Gatherv(local_values_copy.data(), local_nnz, MPI_C_DOUBLE_COMPLEX, res.values.data(), g_counts.data(),
                g_displs.data(), MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Gatherv(local_index_copy.data(), local_nnz, MPI_INT, res.index_col.data(), g_counts.data(), g_displs.data(),
                MPI_INT, 0, MPI_COMM_WORLD);
  } else {
    MPI_Gatherv(res.values.data() + res.index_row[local_start], local_nnz, MPI_C_DOUBLE_COMPLEX, nullptr, nullptr,
                nullptr, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Gatherv(res.index_col.data() + res.index_row[local_start], local_nnz, MPI_INT, nullptr, nullptr, nullptr,
                MPI_INT, 0, MPI_COMM_WORLD);
  }
}

void BroadcastResult(posternak_a_crs_mul_complex_matrix::CRSMatrix &res, int total_nnz) {
  MPI_Bcast(res.values.data(), total_nnz, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(res.index_col.data(), total_nnz, MPI_INT, 0, MPI_COMM_WORLD);
}

}  // namespace

namespace posternak_a_crs_mul_complex_matrix {

PosternakACRSMulComplexMatrixALL::PosternakACRSMulComplexMatrixALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CRSMatrix{};
}

bool PosternakACRSMulComplexMatrixALL::ValidationImpl() {
  const auto &input = GetInput();
  const auto &a = input.first;
  const auto &b = input.second;
  return a.IsValid() && b.IsValid() && a.cols == b.rows;
}

bool PosternakACRSMulComplexMatrixALL::PreProcessingImpl() {
  const auto &input = GetInput();
  const auto &a = input.first;
  const auto &b = input.second;
  auto &res = GetOutput();

  res.rows = a.rows;
  res.cols = b.cols;
  return true;
}

bool PosternakACRSMulComplexMatrixALL::RunImpl() {
  const auto &input = GetInput();
  const auto &a = input.first;
  const auto &b = input.second;
  auto &res = GetOutput();

  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (a.values.empty() || b.values.empty()) {
    res.values.clear();
    res.index_col.clear();
    res.index_row.assign(res.rows + 1, 0);
    return true;
  }

  constexpr double kThreshold = 1e-12;

  int rows_per_proc = res.rows / size;
  int rem = res.rows % size;
  int local_start = (rank * rows_per_proc) + std::min(rank, rem);
  int local_end = local_start + rows_per_proc + (rank < rem ? 1 : 0);
  int local_count = local_end - local_start;

  std::vector<size_t> local_counts(local_count);
  std::vector<size_t> global_counts(res.rows);
  ComputeLocalCounts(a, b, local_counts, local_start, local_count, kThreshold);

  std::vector<int> recv_counts(size);
  std::vector<int> displs(size);
  for (int process = 0; process < size; ++process) {
    int process_start = (process * rows_per_proc) + std::min(process, rem);
    int process_end = process_start + rows_per_proc + (process < rem ? 1 : 0);
    recv_counts[process] = process_end - process_start;
    displs[process] = process_start;
  }

  GatherCountsToRoot(local_counts, local_count, global_counts, recv_counts, displs, rank);

  std::vector<size_t> global_counts_copy = global_counts;
  BuildAndBroadcastStructure(res, global_counts_copy, rank);

  int total_nnz = res.index_row.back();
  res.values.resize(total_nnz);
  res.index_col.resize(total_nnz);

  ComputeLocalRows(a, b, res, local_start, local_count, kThreshold);

  std::vector<int> g_counts(size);
  std::vector<int> g_displs(size);
  PrepareGatherParams(g_counts, g_displs, res, rows_per_proc, rem, size);

  int local_nnz = res.index_row[local_end] - res.index_row[local_start];

  GatherResultData(res, local_start, local_nnz, g_counts, g_displs, rank);
  BroadcastResult(res, total_nnz);

  return res.IsValid();
}

bool PosternakACRSMulComplexMatrixALL::PostProcessingImpl() {
  return GetOutput().IsValid();
}

}  // namespace posternak_a_crs_mul_complex_matrix
