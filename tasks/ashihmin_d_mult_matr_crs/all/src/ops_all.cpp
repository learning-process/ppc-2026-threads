#include "ashihmin_d_mult_matr_crs/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>
#include <tbb/tbb.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <thread>
#include <vector>

#include "ashihmin_d_mult_matr_crs/common/include/common.hpp"
#include "util/include/util.hpp"

namespace ashihmin_d_mult_matr_crs {

AshihminDMultMatrCrsALL::AshihminDMultMatrCrsALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool AshihminDMultMatrCrsALL::ValidationImpl() {
  return GetInput().first.cols == GetInput().second.rows;
}

bool AshihminDMultMatrCrsALL::PreProcessingImpl() {
  auto &matrix_c = GetOutput();
  matrix_c.rows = GetInput().first.rows;
  matrix_c.cols = GetInput().second.cols;
  return true;
}

bool AshihminDMultMatrCrsALL::RunImpl() {
  const auto &matrix_a = GetInput().first;
  const auto &matrix_b = GetInput().second;
  auto &matrix_c = GetOutput();

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int rows_a = matrix_a.rows;

  int base_rows = rows_a / size;
  int rem = rows_a % size;
  int my_start = rank * base_rows + std::min(rank, rem);
  int my_end = my_start + base_rows + (rank < rem ? 1 : 0);
  int my_row_count = my_end - my_start;

  std::vector<std::vector<int>> local_cols(my_row_count);
  std::vector<std::vector<double>> local_vals(my_row_count);

  int thread_count = ppc::util::GetNumThreads();
  std::vector<std::thread> threads;

  auto compute_rows = [&](int start_idx, int end_idx) {
    tbb::parallel_for(start_idx, end_idx, [&](int i) {
      int global_row = my_start + i;
      std::map<int, double> row_accumulator;

      for (int j = matrix_a.row_ptr[global_row]; j < matrix_a.row_ptr[global_row + 1]; ++j) {
        int col_a = matrix_a.col_index[j];
        double val_a = matrix_a.values[j];
        for (int k = matrix_b.row_ptr[col_a]; k < matrix_b.row_ptr[col_a + 1]; ++k) {
          row_accumulator[matrix_b.col_index[k]] += val_a * matrix_b.values[k];
        }
      }

      for (auto it = row_accumulator.begin(); it != row_accumulator.end(); ++it) {
        if (std::abs(it->second) > 1e-15) {
          local_cols[i].push_back(it->first);
          local_vals[i].push_back(it->second);
        }
      }
    });
  };

  int stl_chunk = (my_row_count + thread_count - 1) / thread_count;
  for (int t = 0; t < thread_count; ++t) {
    int s = t * stl_chunk;
    int e = std::min(s + stl_chunk, my_row_count);
    if (s < e) {
      threads.emplace_back(compute_rows, s, e);
    }
  }
  for (auto &th : threads) {
    th.join();
  }

  std::vector<int> my_nnz_per_row(my_row_count);
#pragma omp parallel for
  for (int i = 0; i < my_row_count; ++i) {
    my_nnz_per_row[i] = static_cast<int>(local_cols[i].size());
  }

  std::vector<int> my_flat_cols;
  std::vector<double> my_flat_vals;
  for (int i = 0; i < my_row_count; ++i) {
    my_flat_cols.insert(my_flat_cols.end(), local_cols[i].begin(), local_cols[i].end());
    my_flat_vals.insert(my_flat_vals.end(), local_vals[i].begin(), local_vals[i].end());
  }

  std::vector<int> all_nnz_per_row(rows_a);
  std::vector<int> recv_counts(size);
  std::vector<int> displs(size);

  for (int i = 0; i < size; ++i) {
    recv_counts[i] = (rows_a / size) + (i < (rows_a % size) ? 1 : 0);
    displs[i] = (i == 0) ? 0 : displs[i - 1] + recv_counts[i - 1];
  }

  MPI_Allgatherv(my_nnz_per_row.data(), my_row_count, MPI_INT, all_nnz_per_row.data(), recv_counts.data(),
                 displs.data(), MPI_INT, MPI_COMM_WORLD);

  matrix_c.row_ptr.assign(rows_a + 1, 0);
  for (int i = 0; i < rows_a; ++i) {
    matrix_c.row_ptr[i + 1] = matrix_c.row_ptr[i] + all_nnz_per_row[i];
  }

  int total_nnz = matrix_c.row_ptr[rows_a];
  matrix_c.col_index.resize(total_nnz);
  matrix_c.values.resize(total_nnz);

  std::vector<int> val_recv_counts(size);
  std::vector<int> val_displs(size);
  for (int i = 0; i < size; ++i) {
    val_recv_counts[i] = matrix_c.row_ptr[displs[i] + recv_counts[i]] - matrix_c.row_ptr[displs[i]];
    val_displs[i] = matrix_c.row_ptr[displs[i]];
  }

  MPI_Allgatherv(my_flat_cols.data(), static_cast<int>(my_flat_cols.size()), MPI_INT, matrix_c.col_index.data(),
                 val_recv_counts.data(), val_displs.data(), MPI_INT, MPI_COMM_WORLD);
  MPI_Allgatherv(my_flat_vals.data(), static_cast<int>(my_flat_vals.size()), MPI_DOUBLE, matrix_c.values.data(),
                 val_recv_counts.data(), val_displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);

  return true;
}

bool AshihminDMultMatrCrsALL::PostProcessingImpl() {
  return true;
}

}  // namespace ashihmin_d_mult_matr_crs
