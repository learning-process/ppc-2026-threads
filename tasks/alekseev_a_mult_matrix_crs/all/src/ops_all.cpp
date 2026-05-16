#include "alekseev_a_mult_matrix_crs/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "alekseev_a_mult_matrix_crs/common/include/common.hpp"

#ifdef _MSC_VER
#  define MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
#  define MPI_SIZE_T MPI_UNSIGNED_LONG
#endif

namespace alekseev_a_mult_matrix_crs {

AlekseevAMultMatrixCRSALL::AlekseevAMultMatrixCRSALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool AlekseevAMultMatrixCRSALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int global_valid = 0;

  if (rank == 0) {
    const auto &a = std::get<0>(GetInput());
    const auto &b = std::get<1>(GetInput());
    if (a.cols == b.rows && !a.row_ptr.empty() && !b.row_ptr.empty()) {
      global_valid = 1;
    }
  }
  MPI_Bcast(&global_valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return global_valid == 1;
}

bool AlekseevAMultMatrixCRSALL::PreProcessingImpl() {
  return true;
}

void AlekseevAMultMatrixCRSALL::BroadcastData(CRSMatrix &a, CRSMatrix &b, int rank) {
  std::vector<std::size_t> dims(4);
  if (rank == 0) {
    dims = {a.rows, a.cols, b.rows, b.cols};
  }
  MPI_Bcast(dims.data(), 4, MPI_SIZE_T, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    a.rows = dims[0];
    a.cols = dims[1];
    b.rows = dims[2];
    b.cols = dims[3];
    a.row_ptr.resize(a.rows + 1);
    b.row_ptr.resize(b.rows + 1);
  }

  MPI_Bcast(a.row_ptr.data(), static_cast<int>(a.rows + 1), MPI_SIZE_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(b.row_ptr.data(), static_cast<int>(b.rows + 1), MPI_SIZE_T, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    a.values.resize(a.row_ptr.back());
    a.col_indices.resize(a.row_ptr.back());
    b.values.resize(b.row_ptr.back());
    b.col_indices.resize(b.row_ptr.back());
  }

  MPI_Bcast(a.values.data(), static_cast<int>(a.values.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(a.col_indices.data(), static_cast<int>(a.col_indices.size()), MPI_SIZE_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(b.values.data(), static_cast<int>(b.values.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(b.col_indices.data(), static_cast<int>(b.col_indices.size()), MPI_SIZE_T, 0, MPI_COMM_WORLD);
}

void AlekseevAMultMatrixCRSALL::GatherData(CRSMatrix &out, const InType &input,
                                           const std::vector<std::vector<double>> &local_v,
                                           const std::vector<std::vector<std::size_t>> &local_c,
                                           const std::vector<int> &send_counts, const std::vector<int> &displs,
                                           int rank, int size) {
  int local_n = static_cast<int>(local_v.size());
  std::vector<int> local_row_sizes(local_n);
  std::vector<double> flat_v;
  std::vector<std::size_t> flat_c;

  for (int i = 0; i < local_n; ++i) {
    local_row_sizes[i] = static_cast<int>(local_v[i].size());
    flat_v.insert(flat_v.end(), local_v[i].begin(), local_v[i].end());
    flat_c.insert(flat_c.end(), local_c[i].begin(), local_c[i].end());
  }

  std::vector<int> all_row_sizes;
  if (rank == 0) {
    all_row_sizes.resize(std::get<0>(input).rows);
  }
  MPI_Gatherv(local_row_sizes.data(), local_n, MPI_INT, all_row_sizes.data(), send_counts.data(), displs.data(),
              MPI_INT, 0, MPI_COMM_WORLD);

  int local_total = static_cast<int>(flat_v.size());
  std::vector<int> all_total_counts(size);
  MPI_Gather(&local_total, 1, MPI_INT, all_total_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> v_displs(size, 0);
  std::vector<double> final_v;
  std::vector<std::size_t> final_c;
  int total_elements = 0;

  if (rank == 0) {
    for (int i = 0; i < size; ++i) {
      v_displs[i] = total_elements;
      total_elements += all_total_counts[i];
    }
    final_v.resize(total_elements);
    final_c.resize(total_elements);
  }

  MPI_Gatherv(flat_v.data(), local_total, MPI_DOUBLE, final_v.data(), all_total_counts.data(), v_displs.data(),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(flat_c.data(), local_total, MPI_SIZE_T, final_c.data(), all_total_counts.data(), v_displs.data(),
              MPI_SIZE_T, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    out.rows = std::get<0>(input).rows;
    out.cols = std::get<1>(input).cols;
  }

  MPI_Bcast(&out.rows, 1, MPI_SIZE_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(&out.cols, 1, MPI_SIZE_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(&total_elements, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    final_v.resize(total_elements);
    final_c.resize(total_elements);
    all_row_sizes.resize(out.rows);
  }

  MPI_Bcast(final_v.data(), total_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(final_c.data(), total_elements, MPI_SIZE_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(all_row_sizes.data(), static_cast<int>(out.rows), MPI_INT, 0, MPI_COMM_WORLD);

  out.values = std::move(final_v);
  out.col_indices = std::move(final_c);
  out.row_ptr.assign(out.rows + 1, 0);
  for (std::size_t i = 0; i < out.rows; ++i) {
    out.row_ptr[i + 1] = out.row_ptr[i] + static_cast<std::size_t>(all_row_sizes[i]);
  }
}

bool AlekseevAMultMatrixCRSALL::RunImpl() {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  auto &a = std::get<0>(GetInput());
  auto &b = std::get<1>(GetInput());

  BroadcastData(a, b, rank);

  std::vector<int> send_counts(size);
  std::vector<int> displs(size);
  if (rank == 0) {
    int rows_per_proc = static_cast<int>(a.rows) / size;
    int extra = static_cast<int>(a.rows) % size;
    int offset = 0;
    for (int i = 0; i < size; ++i) {
      send_counts[i] = rows_per_proc + (i < extra ? 1 : 0);
      displs[i] = offset;
      offset += send_counts[i];
    }
  }

  int local_n = 0;
  int local_start = 0;
  MPI_Scatter(send_counts.data(), 1, MPI_INT, &local_n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(displs.data(), 1, MPI_INT, &local_start, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<std::vector<double>> local_v(local_n);
  std::vector<std::vector<std::size_t>> local_c(local_n);

  auto *p_a = &a;
  auto *p_b = &b;
  auto *p_lv = &local_v;
  auto *p_lc = &local_c;

#pragma omp parallel default(none) shared(p_a, p_b, p_lv, p_lc, local_n, local_start)
  {
    std::vector<double> accum(p_b->cols, 0.0);
    std::vector<int> touched_flag(p_b->cols, -1);
    std::vector<std::size_t> touched_cols;
    touched_cols.reserve(p_b->cols);
#pragma omp for schedule(dynamic)
    for (int i = 0; i < local_n; ++i) {
      auto g_row = static_cast<std::size_t>(local_start) + static_cast<std::size_t>(i);
      ProcessRow(g_row, *p_a, *p_b, (*p_lv)[i], (*p_lc)[i], accum, touched_flag, touched_cols);
    }
  }

  GatherData(GetOutput(), GetInput(), local_v, local_c, send_counts, displs, rank, size);
  return true;
}

void AlekseevAMultMatrixCRSALL::ProcessRow(std::size_t i, const CRSMatrix &a, const CRSMatrix &b,
                                           std::vector<double> &temp_v, std::vector<std::size_t> &temp_c,
                                           std::vector<double> &accum, std::vector<int> &touched_flag,
                                           std::vector<std::size_t> &touched_cols) {
  for (std::size_t pos_a = a.row_ptr[i]; pos_a < a.row_ptr[i + 1]; ++pos_a) {
    const std::size_t k = a.col_indices[pos_a];
    const double val_a = a.values[pos_a];

    for (std::size_t pos_b = b.row_ptr[k]; pos_b < b.row_ptr[k + 1]; ++pos_b) {
      const std::size_t j = b.col_indices[pos_b];
      if (std::cmp_not_equal(touched_flag[j], i)) {
        touched_flag[j] = static_cast<int>(i);
        touched_cols.push_back(j);
        accum[j] = 0.0;
      }
      accum[j] += val_a * b.values[pos_b];
    }
  }
  std::ranges::sort(touched_cols);
  for (const auto col : touched_cols) {
    if (std::abs(accum[col]) > 1e-15) {
      temp_v.push_back(accum[col]);
      temp_c.push_back(col);
    }
  }
  touched_cols.clear();
}

bool AlekseevAMultMatrixCRSALL::PostProcessingImpl() {
  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

}  // namespace alekseev_a_mult_matrix_crs
