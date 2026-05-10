#include "../include/ops_mpi.hpp"
#include <mpi.h>
#include <cmath>
#include <vector>
#include <iostream>

namespace nazyrov_a_a_striped_multiplication {

StripedMultiplicationMPI::StripedMultiplicationMPI(const InType &in) : BaseTask() {
  GetInput() = in;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &size_);
}

bool StripedMultiplicationMPI::ValidationImpl() {
  return true;
}

bool StripedMultiplicationMPI::PreProcessingImpl() {
  auto input = GetInput();
  A_ = input.first;
  B_ = input.second;
  n_ = static_cast<int>(std::sqrt(A_.size()));
  
  rows_per_proc_ = n_ / size_;
  if (rank_ < n_ % size_) rows_per_proc_++;
  
  local_C_.resize(rows_per_proc_ * n_, 0.0);
  return true;
}

bool StripedMultiplicationMPI::RunImpl() {
  int start_row = 0;
  for (int i = 0; i < rank_; ++i) {
    start_row += (i < n_ % size_) ? n_ / size_ + 1 : n_ / size_;
  }

  for (int i = 0; i < rows_per_proc_; ++i) {
    int global_row = start_row + i;
    for (int k = 0; k < n_; ++k) {
      double aik = A_[global_row * n_ + k];
      if (aik != 0.0) {
        for (int j = 0; j < n_; ++j) {
          local_C_[i * n_ + j] += aik * B_[k * n_ + j];
        }
      }
    }
  }

  if (rank_ == 0) {
    std::vector<double> global_C(n_ * n_, 0.0);
    for (int i = 0; i < rows_per_proc_; ++i) {
      int global_row = start_row + i;
      for (int j = 0; j < n_; ++j) {
        global_C[global_row * n_ + j] = local_C_[i * n_ + j];
      }
    }
    
    for (int p = 1; p < size_; ++p) {
      int p_rows = (p < n_ % size_) ? n_ / size_ + 1 : n_ / size_;
      std::vector<double> temp(p_rows * n_);
      MPI_Recv(temp.data(), p_rows * n_, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      
      int p_start_row = 0;
      for (int i = 0; i < p; ++i) {
        p_start_row += (i < n_ % size_) ? n_ / size_ + 1 : n_ / size_;
      }
      for (int i = 0; i < p_rows; ++i) {
        for (int j = 0; j < n_; ++j) {
          global_C[(p_start_row + i) * n_ + j] = temp[i * n_ + j];
        }
      }
    }
    GetOutput() = global_C;
  } else {
    MPI_Send(local_C_.data(), rows_per_proc_ * n_, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }

  return true;
}

bool StripedMultiplicationMPI::PostProcessingImpl() {
  return true;
}

}  // namespace nazyrov_a_a_striped_multiplication
