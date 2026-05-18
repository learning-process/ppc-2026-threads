#include "viderman_a_sparse_matrix_mult_crs_complex/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <complex>
#include <vector>

#include "viderman_a_sparse_matrix_mult_crs_complex/common/include/common.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {

VidermanASparseMatrixMultCRSComplexALL::VidermanASparseMatrixMultCRSComplexALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CRSMatrix(0, 0);
}

bool VidermanASparseMatrixMultCRSComplexALL::ValidationImpl() {
  MPI_Initialized(&initialized_pipiline_);

  rank_ = 0;
  if (initialized_pipiline_ != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  }

  int is_valid = 1;
  if (rank_ == 0) {
    const auto &input = GetInput();
    const auto &a = std::get<0>(input);
    const auto &b = std::get<1>(input);

    if (!a.IsValid() || !b.IsValid() || a.cols != b.rows) {
      is_valid = 0;
    }
  }

  if (initialized_pipiline_ != 0) {
    MPI_Bcast(&is_valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }

  return is_valid == 1;
}

bool VidermanASparseMatrixMultCRSComplexALL::PreProcessingImpl() {
  if (initialized_pipiline_ != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size_);
  } else {
    rank_ = 0;
    world_size_ = 1;
  }

  if (rank_ == 0) {
    const auto &input = GetInput();
    const auto &a = std::get<0>(input);
    b_ = std::get<1>(input);

    a_rows_ = a.rows;
    a_cols_ = a.cols;
  }

  if (initialized_pipiline_ != 0) {
    MPI_Bcast(&a_rows_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&a_cols_, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }

  BroadcastData();
  ScatterData();

  return true;
}

bool VidermanASparseMatrixMultCRSComplexALL::RunImpl() {
  MultiplyLocal();

  return true;
}

bool VidermanASparseMatrixMultCRSComplexALL::PostProcessingImpl() {
  GatherData();

  return true;
}

void VidermanASparseMatrixMultCRSComplexALL::MasterScatter() {
  SendDataProcesses();

  int count = a_rows_ / world_size_;
  if (a_rows_ % world_size_ > 0) {
    count++;
  }

  local_a_ = CRSMatrix(count, a_cols_);
  const auto &a = std::get<0>(GetInput());
  if (count > 0 && !a.row_ptr.empty()) {
    int nnz = a.row_ptr[count] - a.row_ptr[0];
    std::copy(a.row_ptr.begin(), a.row_ptr.begin() + count + 1, local_a_.row_ptr.begin());

    local_a_.col_indices.assign(a.col_indices.begin(), a.col_indices.begin() + nnz);
    local_a_.values.assign(a.values.begin(), a.values.begin() + nnz);
  }
}

void VidermanASparseMatrixMultCRSComplexALL::WorkerScatter() {
  int count = 0;
  MPI_Recv(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  local_a_ = CRSMatrix(count, a_cols_);

  if (count > 0) {
    int nnz = 0;
    MPI_Recv(&nnz, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    local_a_.row_ptr.resize(count + 1);
    MPI_Recv(local_a_.row_ptr.data(), count + 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (nnz > 0) {
      local_a_.col_indices.resize(nnz);
      local_a_.values.resize(nnz);

      MPI_Recv(local_a_.col_indices.data(), nnz, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(reinterpret_cast<double *>(local_a_.values.data()), nnz * 2, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }
  }
}

void VidermanASparseMatrixMultCRSComplexALL::ScatterData() {
  if (initialized_pipiline_ == 0 || world_size_ == 1) {
    if (rank_ == 0) {
      local_a_ = std::get<0>(GetInput());
    }

    return;
  }

  if (rank_ == 0) {
    MasterScatter();
  } else {
    WorkerScatter();
  }
}

void VidermanASparseMatrixMultCRSComplexALL::MasterGather() {
  CRSMatrix &out = GetOutput();
  out = CRSMatrix(a_rows_, b_.cols);

  int rows_to_copy = std::min(local_c_.rows, a_rows_);
  for (int i = 0; i <= rows_to_copy; ++i) {
    out.row_ptr[i] = local_c_.row_ptr[i];
  }

  out.col_indices = local_c_.col_indices;
  out.values = local_c_.values;

  ReceiveDataProcesses(out);
}

void VidermanASparseMatrixMultCRSComplexALL::WorkerGather() {
  int send_rows = local_c_.rows;
  int send_nnz = static_cast<int>(local_c_.values.size());

  MPI_Send(&send_rows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

  if (send_rows > 0) {
    MPI_Send(&send_nnz, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
    MPI_Send(local_c_.row_ptr.data(), send_rows + 1, MPI_INT, 0, 2, MPI_COMM_WORLD);

    if (send_nnz > 0) {
      MPI_Send(local_c_.col_indices.data(), send_nnz, MPI_INT, 0, 3, MPI_COMM_WORLD);
      MPI_Send(reinterpret_cast<const double *>(local_c_.values.data()), send_nnz * 2, MPI_DOUBLE, 0, 4,
               MPI_COMM_WORLD);
    }
  }
}

void VidermanASparseMatrixMultCRSComplexALL::GatherData() {
  if (initialized_pipiline_ == 0 || world_size_ == 1) {
    if (rank_ == 0) {
      GetOutput() = local_c_;
    }

    return;
  }

  if (rank_ == 0) {
    MasterGather();
  } else {
    WorkerGather();
  }

  CRSMatrix &out = GetOutput();

  std::vector<int> g_info(3, 0);
  if (rank_ == 0) {
    g_info[0] = out.rows;
    g_info[1] = out.cols;
    g_info[2] = static_cast<int>(out.values.size());
  }

  MPI_Bcast(g_info.data(), 3, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank_ != 0) {
    out.rows = g_info[0];
    out.cols = g_info[1];
    out.row_ptr.assign(out.rows + 1, 0);
    out.col_indices.resize(g_info[2]);
    out.values.resize(g_info[2]);
  }

  if (out.rows >= 0) {
    MPI_Bcast(out.row_ptr.data(), out.rows + 1, MPI_INT, 0, MPI_COMM_WORLD);
  }

  if (g_info[2] > 0) {
    MPI_Bcast(out.col_indices.data(), g_info[2], MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(reinterpret_cast<double *>(out.values.data()), g_info[2] * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
}

void VidermanASparseMatrixMultCRSComplexALL::SendDataProcesses() {
  const auto &a = std::get<0>(GetInput());
  int rows_per_proc = a_rows_ / world_size_;
  int remainder = a_rows_ % world_size_;

  for (int process = 1; process < world_size_; ++process) {
    int start = (process * rows_per_proc) + std::min(process, remainder);
    int count = rows_per_proc + (process < remainder ? 1 : 0);
    MPI_Send(&count, 1, MPI_INT, process, 0, MPI_COMM_WORLD);

    if (count > 0) {
      int nnz = a.row_ptr[start + count] - a.row_ptr[start];
      MPI_Send(&nnz, 1, MPI_INT, process, 1, MPI_COMM_WORLD);

      std::vector<int> p_row_ptr(count + 1);

      int offset = a.row_ptr[start];
      for (int i = 0; i <= count; ++i) {
        p_row_ptr[i] = a.row_ptr[start + i] - offset;
      }

      MPI_Send(p_row_ptr.data(), count + 1, MPI_INT, process, 2, MPI_COMM_WORLD);

      if (nnz > 0) {
        MPI_Send(&a.col_indices[offset], nnz, MPI_INT, process, 3, MPI_COMM_WORLD);
        MPI_Send(reinterpret_cast<const double *>(&a.values[offset]), nnz * 2, MPI_DOUBLE, process, 4, MPI_COMM_WORLD);
      }
    }
  }
}

void VidermanASparseMatrixMultCRSComplexALL::ReceiveDataProcesses(CRSMatrix &out) const {
  int current_global_row = local_c_.rows;
  for (int process = 1; process < world_size_; ++process) {
    int p_rows = 0;
    MPI_Recv(&p_rows, 1, MPI_INT, process, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (p_rows > 0) {
      int p_nnz = 0;
      MPI_Recv(&p_nnz, 1, MPI_INT, process, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      std::vector<int> p_ptr(p_rows + 1);
      MPI_Recv(p_ptr.data(), p_rows + 1, MPI_INT, process, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      int offset = static_cast<int>(out.values.size());

      if (p_nnz > 0) {
        std::vector<int> p_cols(p_nnz);
        std::vector<Complex> p_vals(p_nnz);

        MPI_Recv(p_cols.data(), p_nnz, MPI_INT, process, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(reinterpret_cast<double *>(p_vals.data()), p_nnz * 2, MPI_DOUBLE, process, 4, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        out.col_indices.insert(out.col_indices.end(), p_cols.begin(), p_cols.end());
        out.values.insert(out.values.end(), p_vals.begin(), p_vals.end());
      }
      for (int i = 1; i <= p_rows; ++i) {
        if (current_global_row + i <= a_rows_) {
          out.row_ptr[current_global_row + i] = offset + p_ptr[i];
        }
      }

      current_global_row += p_rows;
    }
  }
}

void VidermanASparseMatrixMultCRSComplexALL::BroadcastData() {
  if (initialized_pipiline_ == 0) {
    return;
  }

  MPI_Bcast(&b_.rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&b_.cols, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int nnz = (rank_ == 0) ? static_cast<int>(b_.values.size()) : 0;
  MPI_Bcast(&nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank_ != 0) {
    b_.row_ptr.assign(b_.rows + 1, 0);
    b_.col_indices.resize(nnz);
    b_.values.resize(nnz);
  }

  if (!b_.row_ptr.empty()) {
    MPI_Bcast(b_.row_ptr.data(), b_.rows + 1, MPI_INT, 0, MPI_COMM_WORLD);
  }

  if (nnz > 0) {
    MPI_Bcast(b_.col_indices.data(), nnz, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(reinterpret_cast<double *>(b_.values.data()), nnz * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
}

void VidermanASparseMatrixMultCRSComplexALL::ComputeMultiply(int i, std::vector<int> &cols,
                                                             std::vector<Complex> &vals) {
  std::vector<Complex> accumulator(b_.cols, Complex(0, 0));
  bool has_data = false;

  for (int j = local_a_.row_ptr[i]; j < local_a_.row_ptr[i + 1]; ++j) {
    int a_col = local_a_.col_indices[j];
    Complex a_val = local_a_.values[j];

    for (int k = b_.row_ptr[a_col]; k < b_.row_ptr[a_col + 1]; ++k) {
      accumulator[b_.col_indices[k]] += a_val * b_.values[k];

      has_data = true;
    }
  }
  if (has_data) {
    for (int j = 0; j < b_.cols; ++j) {
      if (std::norm(accumulator[j]) > 1e-18) {
        cols.push_back(j);
        vals.push_back(accumulator[j]);
      }
    }
  }
}

void VidermanASparseMatrixMultCRSComplexALL::MultiplyLocal() {
  local_c_ = CRSMatrix(local_a_.rows, b_.cols);
  if (local_a_.rows <= 0 || b_.cols <= 0) {
    return;
  }

  std::vector<std::vector<int>> row_cols(local_a_.rows);
  std::vector<std::vector<Complex>> row_vals(local_a_.rows);

#pragma omp parallel for default(none) shared(row_cols, row_vals)
  for (int i = 0; i < local_a_.rows; ++i) {
    ComputeMultiply(i, row_cols[i], row_vals[i]);
  }

  for (int i = 0; i < local_a_.rows; ++i) {
    local_c_.col_indices.insert(local_c_.col_indices.end(), row_cols[i].begin(), row_cols[i].end());
    local_c_.values.insert(local_c_.values.end(), row_vals[i].begin(), row_vals[i].end());
    local_c_.row_ptr[i + 1] = static_cast<int>(local_c_.values.size());
  }
}

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
