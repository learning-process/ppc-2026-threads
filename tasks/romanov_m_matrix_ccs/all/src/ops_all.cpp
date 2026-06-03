#include "romanov_m_matrix_ccs/all/include/ops_all.hpp"

#include <mpi.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "romanov_m_matrix_ccs/common/include/common.hpp"

namespace romanov_m_matrix_ccs {

namespace {

void BroadcastSizeTVector(std::vector<size_t> &data, size_t size, int rank) {
  std::vector<uint64_t> buffer;
  if (rank == 0) {
    buffer.assign(data.begin(), data.end());
  } else {
    buffer.resize(size);
  }

  MPI_Bcast(static_cast<void *>(buffer.data()), static_cast<int>(size), MPI_UINT64_T, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    data.assign(buffer.begin(), buffer.end());
  }
}

void SendSizeTVector(const std::vector<size_t> &data) {
  std::vector<uint64_t> buffer(data.begin(), data.end());
  MPI_Send(static_cast<const void *>(buffer.data()), static_cast<int>(buffer.size()), MPI_UINT64_T, 0, 2,
           MPI_COMM_WORLD);
}

void RecvSizeTVector(std::vector<size_t> &data, int nnz, int proc) {
  std::vector<uint64_t> buffer(static_cast<size_t>(nnz));
  MPI_Recv(static_cast<void *>(buffer.data()), nnz, MPI_UINT64_T, proc, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  data.assign(buffer.begin(), buffer.end());
}

}  // namespace

RomanovMMatrixCCSALL::RomanovMMatrixCCSALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool RomanovMMatrixCCSALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int res = 0;
  if (rank == 0) {
    const auto &left = GetInput().first;
    const auto &right = GetInput().second;
    res = (left.cols_num == right.rows_num && left.cols_num > 0) ? 1 : 0;
  }
  MPI_Bcast(&res, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return res == 1;
}

bool RomanovMMatrixCCSALL::PreProcessingImpl() {
  return true;
}

void RomanovMMatrixCCSALL::MultiplyColumn(size_t col_index, const MatrixCCS &a, const MatrixCCS &b,
                                          std::vector<double> &temp_v, std::vector<size_t> &temp_r) {
  std::vector<double> accumulator(a.rows_num, 0.0);
  std::vector<bool> row_mask(a.rows_num, false);
  std::vector<size_t> active_rows;

  for (size_t kb = b.col_ptrs[col_index]; kb < b.col_ptrs[col_index + 1]; ++kb) {
    size_t k = b.row_inds[kb];
    double v_b = b.vals[kb];
    for (size_t ka = a.col_ptrs[k]; ka < a.col_ptrs[k + 1]; ++ka) {
      size_t r_idx = a.row_inds[ka];
      if (!row_mask[r_idx]) {
        row_mask[r_idx] = true;
        active_rows.push_back(r_idx);
      }
      accumulator[r_idx] += a.vals[ka] * v_b;
    }
  }

  std::ranges::sort(active_rows);
  for (size_t r_idx : active_rows) {
    if (std::abs(accumulator[r_idx]) > 1e-12) {
      temp_v.push_back(accumulator[r_idx]);
      temp_r.push_back(r_idx);
    }
  }
}

void RomanovMMatrixCCSALL::SyncMatrixData(int rank, MatrixCCS &a, MatrixCCS &b) {
  std::array<uint64_t, 3> dims{};
  if (rank == 0) {
    dims[0] = static_cast<uint64_t>(a.rows_num);
    dims[1] = static_cast<uint64_t>(a.cols_num);
    dims[2] = static_cast<uint64_t>(b.cols_num);
  }
  MPI_Bcast(static_cast<void *>(dims.data()), 3, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    a.rows_num = static_cast<size_t>(dims[0]);
    a.cols_num = static_cast<size_t>(dims[1]);
    b.rows_num = static_cast<size_t>(dims[1]);
    b.cols_num = static_cast<size_t>(dims[2]);
    a.col_ptrs.resize(a.cols_num + 1);
    b.col_ptrs.resize(b.cols_num + 1);
  }

  BroadcastSizeTVector(a.col_ptrs, a.cols_num + 1, rank);
  size_t a_nnz = (rank == 0) ? a.vals.size() : a.col_ptrs[a.cols_num];
  if (rank != 0) {
    a.row_inds.resize(a_nnz);
    a.vals.resize(a_nnz);
  }
  BroadcastSizeTVector(a.row_inds, a_nnz, rank);
  if (a_nnz > 0) {
    MPI_Bcast(static_cast<void *>(a.vals.data()), static_cast<int>(a_nnz), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  BroadcastSizeTVector(b.col_ptrs, b.cols_num + 1, rank);
  size_t b_nnz = (rank == 0) ? b.vals.size() : b.col_ptrs[b.cols_num];
  if (rank != 0) {
    b.row_inds.resize(b_nnz);
    b.vals.resize(b_nnz);
  }
  BroadcastSizeTVector(b.row_inds, b_nnz, rank);
  if (b_nnz > 0) {
    MPI_Bcast(static_cast<void *>(b.vals.data()), static_cast<int>(b_nnz), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
}

void RomanovMMatrixCCSALL::MasterCollect(int size, int chunk, int remainder, std::vector<std::vector<double>> &all_v,
                                         std::vector<std::vector<size_t>> &all_r) {
  for (int proc_idx = 1; proc_idx < size; ++proc_idx) {
    int p_start = (proc_idx * chunk) + std::min(proc_idx, remainder);
    int p_cols = chunk + ((proc_idx < remainder) ? 1 : 0);
    for (int i = 0; i < p_cols; ++i) {
      int nnz = 0;
      MPI_Recv(&nnz, 1, MPI_INT, proc_idx, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      auto target_idx = static_cast<size_t>(p_start) + i;
      all_v[target_idx].resize(static_cast<size_t>(nnz));
      all_r[target_idx].resize(static_cast<size_t>(nnz));
      if (nnz > 0) {
        MPI_Recv(static_cast<void *>(all_v[target_idx].data()), nnz, MPI_DOUBLE, proc_idx, 1, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        RecvSizeTVector(all_r[target_idx], nnz, proc_idx);
      }
    }
  }
}

void RomanovMMatrixCCSALL::WorkerSend(int local_count, std::vector<std::vector<double>> &local_v,
                                      std::vector<std::vector<size_t>> &local_r) {
  for (int i = 0; i < local_count; ++i) {
    auto idx = static_cast<size_t>(i);
    int nnz = static_cast<int>(local_v[idx].size());
    MPI_Send(&nnz, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    if (nnz > 0) {
      MPI_Send(static_cast<const void *>(local_v[idx].data()), nnz, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
      SendSizeTVector(local_r[idx]);
    }
  }
}

bool RomanovMMatrixCCSALL::RunImpl() {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MatrixCCS a_mat = (rank == 0) ? GetInput().first : MatrixCCS();
  MatrixCCS b_mat = (rank == 0) ? GetInput().second : MatrixCCS();
  SyncMatrixData(rank, a_mat, b_mat);

  int total_cols = static_cast<int>(b_mat.cols_num);
  int chunk = total_cols / size;
  int remainder = total_cols % size;
  int start_col = (rank * chunk) + std::min(rank, remainder);
  int end_col = start_col + (chunk + ((rank < remainder) ? 1 : 0));
  int local_count = std::max(0, end_col - start_col);

  std::vector<std::vector<double>> local_v(static_cast<size_t>(local_count));
  std::vector<std::vector<size_t>> local_r(static_cast<size_t>(local_count));

  if (local_count > 0) {
    tbb::parallel_for(0, local_count, [&](int i) {
      MultiplyColumn(static_cast<size_t>(start_col) + i, a_mat, b_mat, local_v[static_cast<size_t>(i)],
                     local_r[static_cast<size_t>(i)]);
    });
  }

  auto &c_mat = GetOutput();
  if (rank == 0) {
    std::vector<std::vector<double>> all_v(static_cast<size_t>(total_cols));
    std::vector<std::vector<size_t>> all_r(static_cast<size_t>(total_cols));
    for (int i = 0; i < local_count; ++i) {
      auto global_idx = static_cast<size_t>(start_col) + i;
      auto l_idx = static_cast<size_t>(i);
      all_v[global_idx] = std::move(local_v[l_idx]);
      all_r[global_idx] = std::move(local_r[l_idx]);
    }
    MasterCollect(size, chunk, remainder, all_v, all_r);

    c_mat.rows_num = a_mat.rows_num;
    c_mat.cols_num = static_cast<size_t>(total_cols);
    c_mat.col_ptrs.assign(c_mat.cols_num + 1, 0);
    for (size_t j = 0; j < c_mat.cols_num; ++j) {
      c_mat.col_ptrs[j + 1] = c_mat.col_ptrs[j] + all_v[j].size();
      c_mat.vals.insert(c_mat.vals.end(), all_v[j].begin(), all_v[j].end());
      c_mat.row_inds.insert(c_mat.row_inds.end(), all_r[j].begin(), all_r[j].end());
    }
    c_mat.nnz = c_mat.vals.size();
  } else {
    WorkerSend(local_count, local_v, local_r);
  }

  std::array<uint64_t, 3> final_dims{};
  if (rank == 0) {
    final_dims[0] = static_cast<uint64_t>(c_mat.rows_num);
    final_dims[1] = static_cast<uint64_t>(c_mat.cols_num);
    final_dims[2] = static_cast<uint64_t>(c_mat.nnz);
  }
  MPI_Bcast(static_cast<void *>(final_dims.data()), 3, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    c_mat.rows_num = static_cast<size_t>(final_dims[0]);
    c_mat.cols_num = static_cast<size_t>(final_dims[1]);
    c_mat.nnz = static_cast<size_t>(final_dims[2]);
    c_mat.col_ptrs.resize(c_mat.cols_num + 1);
    c_mat.row_inds.resize(c_mat.nnz);
    c_mat.vals.resize(c_mat.nnz);
  }

  BroadcastSizeTVector(c_mat.col_ptrs, c_mat.cols_num + 1, rank);
  BroadcastSizeTVector(c_mat.row_inds, c_mat.nnz, rank);
  if (c_mat.nnz > 0) {
    MPI_Bcast(static_cast<void *>(c_mat.vals.data()), static_cast<int>(c_mat.nnz), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

bool RomanovMMatrixCCSALL::PostProcessingImpl() {
  return true;
}

}  // namespace romanov_m_matrix_ccs
