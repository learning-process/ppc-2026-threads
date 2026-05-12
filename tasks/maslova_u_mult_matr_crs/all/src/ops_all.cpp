#include "maslova_u_mult_matr_crs/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include "maslova_u_mult_matr_crs/common/include/common.hpp"
#include "util/include/util.hpp"

namespace maslova_u_mult_matr_crs {

void MaslovaUMultMatrALL::SortVector(std::vector<int> &vec) {  // для прохождения clang-tidy
  std::ranges::sort(vec);
}

MaslovaUMultMatrALL::MaslovaUMultMatrALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    GetInput() = in;
  }
}

bool MaslovaUMultMatrALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int ok = 0;
  if (rank == 0) {
    const auto &input = GetInput();
    const auto &matrix_a = std::get<0>(input);
    const auto &matrix_b = std::get<1>(input);
    if (matrix_a.cols == matrix_b.rows && matrix_a.rows > 0 && matrix_b.cols > 0) {
      ok = 1;
    }
  }
  MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return ok == 1;
}

bool MaslovaUMultMatrALL::PreProcessingImpl() {
  return true;
}

void MaslovaUMultMatrALL::BroadcastCRSMatrix(CRSMatrix &m, int root, int r, int c) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::array<int, 2> counts = {0, 0};
  if (rank == root) {
    counts[0] = static_cast<int>(m.row_ptr.size());
    counts[1] = static_cast<int>(m.values.size());
  }
  MPI_Bcast(counts.data(), 2, MPI_INT, root, MPI_COMM_WORLD);
  if (rank != root) {
    m.row_ptr.resize(counts[0]);
    m.values.resize(counts[1]);
    m.col_ind.resize(counts[1]);
  }
  if (counts[0] > 0) {
    MPI_Bcast(m.row_ptr.data(), counts[0], MPI_INT, root, MPI_COMM_WORLD);
  }
  if (counts[1] > 0) {
    MPI_Bcast(m.values.data(), counts[1], MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(m.col_ind.data(), counts[1], MPI_INT, root, MPI_COMM_WORLD);
  }
  m.rows = r;
  m.cols = c;
}

void MaslovaUMultMatrALL::ComputeLocalPart(const CRSMatrix &a, const CRSMatrix &b, int start_row, int local_rows,
                                           std::vector<int> &local_nnz, std::vector<double> &flat_values,
                                           std::vector<int> &flat_cols) {
  if (local_rows <= 0) {
    return;
  }
  std::vector<std::vector<double>> t_vals(local_rows);
  std::vector<std::vector<int>> t_cols(local_rows);

#pragma omp parallel default(none) shared(a, b, start_row, local_rows, local_nnz, t_vals, t_cols) \
    num_threads(ppc::util::GetNumThreads())
  {
    std::vector<int> marker(b.cols, -1);
    std::vector<double> acc(b.cols, 0.0);
    std::vector<int> used;
#pragma omp for schedule(dynamic)
    for (int i = 0; i < local_rows; ++i) {
      int g_row = start_row + i;
      used.clear();
      for (int j = a.row_ptr[g_row]; j < a.row_ptr[g_row + 1]; ++j) {
        int col_a = a.col_ind[j];
        double val_a = a.values[j];
        for (int k = b.row_ptr[col_a]; k < b.row_ptr[col_a + 1]; ++k) {
          int col_b = b.col_ind[k];
          if (marker[col_b] != i) {
            marker[col_b] = i;
            used.push_back(col_b);
            acc[col_b] = val_a * b.values[k];
          } else {
            acc[col_b] += val_a * b.values[k];
          }
        }
      }
      local_nnz[i] = static_cast<int>(used.size());
      SortVector(used);
      for (int col : used) {
        t_vals[i].push_back(acc[col]);
        t_cols[i].push_back(col);
        acc[col] = 0.0;
      }
    }
  }
  for (int i = 0; i < local_rows; ++i) {
    flat_values.insert(flat_values.end(), t_vals[i].begin(), t_vals[i].end());
    flat_cols.insert(flat_cols.end(), t_cols[i].begin(), t_cols[i].end());
  }
}

void MaslovaUMultMatrALL::GatherResults(int rank, int size, int a_rows, int b_cols, int local_rows, CRSMatrix &c,
                                        const std::vector<int> &local_nnz, const std::vector<double> &flat_values,
                                        const std::vector<int> &flat_cols) {
  int local_nnz_total = static_cast<int>(flat_values.size());
  std::vector<int> all_nnz_counts(size);
  std::vector<int> all_row_counts(size);

  MPI_Gather(&local_nnz_total, 1, MPI_INT, all_nnz_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&local_rows, 1, MPI_INT, all_row_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    c.rows = a_rows;
    c.cols = b_cols;
    std::vector<int> nnz_displs(size, 0);
    std::vector<int> row_offsets(size, 0);
    int total_nnz = 0;
    for (int i = 0; i < size; ++i) {
      nnz_displs[i] = total_nnz;
      total_nnz += all_nnz_counts[i];
      if (i > 0) {
        row_offsets[i] = row_offsets[i - 1] + all_row_counts[i - 1];
      }
    }
    c.values.resize(total_nnz);
    c.col_ind.resize(total_nnz);
    c.row_ptr.assign(static_cast<size_t>(a_rows) + 1, 0);

    MPI_Gatherv(flat_values.data(), local_nnz_total, MPI_DOUBLE, c.values.data(), all_nnz_counts.data(),
                nnz_displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(flat_cols.data(), local_nnz_total, MPI_INT, c.col_ind.data(), all_nnz_counts.data(), nnz_displs.data(),
                MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> all_nnz_per_row(a_rows);
    MPI_Gatherv(local_nnz.data(), local_rows, MPI_INT, all_nnz_per_row.data(), all_row_counts.data(),
                row_offsets.data(), MPI_INT, 0, MPI_COMM_WORLD);
    for (int i = 0; i < a_rows; ++i) {
      c.row_ptr[i + 1] = c.row_ptr[i] + all_nnz_per_row[i];
    }
  } else {
    MPI_Gatherv(flat_values.data(), local_nnz_total, MPI_DOUBLE, nullptr, nullptr, nullptr, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
    MPI_Gatherv(flat_cols.data(), local_nnz_total, MPI_INT, nullptr, nullptr, nullptr, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(local_nnz.data(), local_rows, MPI_INT, nullptr, nullptr, nullptr, MPI_INT, 0, MPI_COMM_WORLD);
  }
  BroadcastCRSMatrix(c, 0, a_rows, b_cols);
}

bool MaslovaUMultMatrALL::RunImpl() {
  int size = 0;
  int rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::array<int, 3> dims = {0, 0, 0};
  if (rank == 0) {
    const auto &a_in = std::get<0>(GetInput());
    const auto &b_in = std::get<1>(GetInput());
    dims[0] = a_in.rows;
    dims[1] = a_in.cols;
    dims[2] = b_in.cols;
  }
  MPI_Bcast(dims.data(), 3, MPI_INT, 0, MPI_COMM_WORLD);

  CRSMatrix a;
  CRSMatrix b;
  if (rank == 0) {
    a = std::get<0>(GetInput());
    b = std::get<1>(GetInput());
  }
  BroadcastCRSMatrix(a, 0, dims[0], dims[1]);
  BroadcastCRSMatrix(b, 0, dims[1], dims[2]);

  int part = dims[0] / size;
  int rem = dims[0] % size;
  int start_row = (rank * part) + std::min(rank, rem);
  int local_rows = part + (rank < rem ? 1 : 0);

  std::vector<int> local_nnz(local_rows);
  std::vector<double> flat_values;
  std::vector<int> flat_cols;

  ComputeLocalPart(a, b, start_row, local_rows, local_nnz, flat_values, flat_cols);
  GatherResults(rank, size, dims[0], dims[2], local_rows, GetOutput(), local_nnz, flat_values, flat_cols);

  return true;
}

bool MaslovaUMultMatrALL::PostProcessingImpl() {
  return true;
}

}  // namespace maslova_u_mult_matr_crs
