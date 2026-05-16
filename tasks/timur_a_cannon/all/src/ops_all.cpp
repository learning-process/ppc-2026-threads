#include "timur_a_cannon/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <utility>
#include <vector>

namespace timur_a_cannon {

namespace {

using Matrix = std::vector<std::vector<double>>;

void CopyBlocksForStep(const Matrix &src_a, const Matrix &src_b, int b_size, int global_i, int shift, int j,
                       Matrix &block_a, Matrix &block_b) {
  for (int row = 0; row < b_size; ++row) {
    for (int col = 0; col < b_size; ++col) {
      block_a[row][col] = src_a[(global_i * b_size) + row][(shift * b_size) + col];
      block_b[row][col] = src_b[(shift * b_size) + row][(j * b_size) + col];
    }
  }
}

void ScatterBlockIntoResult(Matrix &local_result, const Matrix &block_c, int local_i, int j, int b_size) {
  for (int row = 0; row < b_size; ++row) {
    for (int col = 0; col < b_size; ++col) {
      local_result[(local_i * b_size) + row][(j * b_size) + col] = block_c[row][col];
    }
  }
}

std::vector<double> FlattenMatrix(const Matrix &matrix) {
  const std::size_t rows = matrix.size();
  const std::size_t cols = rows == 0 ? 0 : matrix[0].size();
  std::vector<double> flat(rows * cols);

  for (std::size_t row = 0; row < rows; ++row) {
    std::copy(matrix[row].begin(), matrix[row].end(), flat.begin() + static_cast<std::ptrdiff_t>(row * cols));
  }

  return flat;
}

Matrix UnflattenMatrix(const std::vector<double> &flat, std::size_t rows, std::size_t cols) {
  Matrix matrix(rows, std::vector<double>(cols));

  for (std::size_t row = 0; row < rows; ++row) {
    const std::ptrdiff_t begin_idx = (static_cast<std::ptrdiff_t>(row) * static_cast<std::ptrdiff_t>(cols));
    const std::ptrdiff_t end_idx = (static_cast<std::ptrdiff_t>(row + 1) * static_cast<std::ptrdiff_t>(cols));
    std::copy(flat.begin() + begin_idx, flat.begin() + end_idx, matrix[row].begin());
  }

  return matrix;
}

std::pair<std::vector<int>, std::vector<int>> BuildGatherLayout(int size, int base_block_rows, int extra_block_rows,
                                                                int b_size, int n) {
  std::vector<int> recv_counts(size);
  std::vector<int> displs(size);
  int offset = 0;
  for (int proc = 0; proc < size; ++proc) {
    const int proc_block_rows = base_block_rows + (proc < extra_block_rows ? 1 : 0);
    recv_counts[proc] = proc_block_rows * b_size * n;
    displs[proc] = offset;
    offset += recv_counts[proc];
  }
  return {recv_counts, displs};
}

}  // namespace

TimurACannonMatrixMultiplicationALL::TimurACannonMatrixMultiplicationALL(
    const std::tuple<int, std::vector<std::vector<double>>, std::vector<std::vector<double>>> &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TimurACannonMatrixMultiplicationALL::ValidationImpl() {
  const auto &input = GetInput();
  const int b_size = std::get<0>(input);
  const auto &mat_a = std::get<1>(input);
  const auto &mat_b = std::get<2>(input);

  if (b_size <= 0 || mat_a.empty() || mat_b.empty()) {
    return false;
  }

  const std::size_t n = mat_a.size();
  if (mat_b.size() != n || (n % static_cast<std::size_t>(b_size) != 0)) {
    return false;
  }

  const auto is_square_n = [n](const Matrix &matrix) {
    return std::ranges::all_of(matrix, [n](const std::vector<double> &row) { return row.size() == n; });
  };

  return is_square_n(mat_a) && is_square_n(mat_b);
}

bool TimurACannonMatrixMultiplicationALL::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

void TimurACannonMatrixMultiplicationALL::BlockMultiplyAccumulate(const std::vector<std::vector<double>> &a,
                                                                  const std::vector<std::vector<double>> &b,
                                                                  std::vector<std::vector<double>> &c, int b_size) {
  for (int i = 0; i < b_size; ++i) {
    for (int k = 0; k < b_size; ++k) {
      const double temp = a[i][k];
      for (int j = 0; j < b_size; ++j) {
        c[i][j] += temp * b[k][j];
      }
    }
  }
}

std::vector<std::vector<double>> TimurACannonMatrixMultiplicationALL::ComputeLocalResult(const Matrix &src_a,
                                                                                         const Matrix &src_b,
                                                                                         int b_size, int grid_sz,
                                                                                         int block_row_start,
                                                                                         int local_block_rows, int n) {
  Matrix local_result(static_cast<std::size_t>(local_block_rows) * static_cast<std::size_t>(b_size),
                      std::vector<double>(static_cast<std::size_t>(n), 0.0));

#pragma omp parallel for default(none) \
    shared(local_result, src_a, src_b, b_size, grid_sz, block_row_start, local_block_rows)
  for (int local_i = 0; local_i < local_block_rows; ++local_i) {
    for (int j = 0; j < grid_sz; ++j) {
      Matrix block_c(b_size, std::vector<double>(b_size, 0.0));
      const int global_i = block_row_start + local_i;

      for (int step = 0; step < grid_sz; ++step) {
        const int shift = (global_i + j + step) % grid_sz;
        Matrix block_a(b_size, std::vector<double>(b_size));
        Matrix block_b(b_size, std::vector<double>(b_size));
        CopyBlocksForStep(src_a, src_b, b_size, global_i, shift, j, block_a, block_b);
        BlockMultiplyAccumulate(block_a, block_b, block_c, b_size);
      }

      ScatterBlockIntoResult(local_result, block_c, local_i, j, b_size);
    }
  }

  return local_result;
}

bool TimurACannonMatrixMultiplicationALL::RunImpl() {
  const auto &input = GetInput();
  const int b_size = std::get<0>(input);
  Matrix src_a = std::get<1>(input);
  Matrix src_b = std::get<2>(input);
  const int n = static_cast<int>(src_a.size());
  const int grid_sz = n / b_size;
  const int total_elems = n * n;

  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<double> flat_a = FlattenMatrix(src_a);
  std::vector<double> flat_b = FlattenMatrix(src_b);

  MPI_Bcast(flat_a.data(), total_elems, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(flat_b.data(), total_elems, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  src_a = UnflattenMatrix(flat_a, static_cast<std::size_t>(n), static_cast<std::size_t>(n));
  src_b = UnflattenMatrix(flat_b, static_cast<std::size_t>(n), static_cast<std::size_t>(n));

  const int base_block_rows = grid_sz / size;
  const int extra_block_rows = grid_sz % size;
  const int local_block_rows = base_block_rows + (rank < extra_block_rows ? 1 : 0);
  const int block_row_start = (rank * base_block_rows) + std::min(rank, extra_block_rows);

  Matrix local_result = ComputeLocalResult(src_a, src_b, b_size, grid_sz, block_row_start, local_block_rows, n);

  std::vector<double> local_flat = FlattenMatrix(local_result);
  auto [recv_counts, displs] = BuildGatherLayout(size, base_block_rows, extra_block_rows, b_size, n);

  std::vector<double> global_flat(total_elems);
  MPI_Allgatherv(local_flat.data(), static_cast<int>(local_flat.size()), MPI_DOUBLE, global_flat.data(),
                 recv_counts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);

  GetOutput() = UnflattenMatrix(global_flat, static_cast<std::size_t>(n), static_cast<std::size_t>(n));
  return true;
}

bool TimurACannonMatrixMultiplicationALL::PostProcessingImpl() {
  return true;
}

}  // namespace timur_a_cannon
