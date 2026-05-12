#include "timur_a_cannon/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <vector>

namespace timur_a_cannon {

namespace {

using Matrix = std::vector<std::vector<double>>;

std::vector<double> FlattenMatrix(const Matrix &matrix) {
  const std::size_t rows = matrix.size();
  const std::size_t cols = rows == 0 ? 0 : matrix[0].size();
  std::vector<double> flat(rows * cols);

  for (std::size_t row = 0; row < rows; ++row) {
    std::copy(matrix[row].begin(), matrix[row].end(), flat.begin() + static_cast<std::ptrdiff_t>(row * cols));
  }

  return flat;
}

Matrix UnflattenMatrix(const std::vector<double> &flat, int rows, int cols) {
  Matrix matrix(rows, std::vector<double>(cols));

  for (int row = 0; row < rows; ++row) {
    std::copy(flat.begin() + static_cast<std::ptrdiff_t>(row * cols),
              flat.begin() + static_cast<std::ptrdiff_t>((row + 1) * cols), matrix[row].begin());
  }

  return matrix;
}

}  // namespace

TimurACannonMatrixMultiplicationALL::TimurACannonMatrixMultiplicationALL(const InType &in) {
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

  src_a = UnflattenMatrix(flat_a, n, n);
  src_b = UnflattenMatrix(flat_b, n, n);

  const int base_block_rows = grid_sz / size;
  const int extra_block_rows = grid_sz % size;
  const int local_block_rows = base_block_rows + (rank < extra_block_rows ? 1 : 0);
  const int block_row_start = (rank * base_block_rows) + std::min(rank, extra_block_rows);

  Matrix local_result(local_block_rows * b_size, std::vector<double>(n, 0.0));

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

        for (int row = 0; row < b_size; ++row) {
          for (int col = 0; col < b_size; ++col) {
            block_a[row][col] = src_a[(global_i * b_size) + row][(shift * b_size) + col];
            block_b[row][col] = src_b[(shift * b_size) + row][(j * b_size) + col];
          }
        }

        BlockMultiplyAccumulate(block_a, block_b, block_c, b_size);
      }

      for (int row = 0; row < b_size; ++row) {
        for (int col = 0; col < b_size; ++col) {
          local_result[(local_i * b_size) + row][(j * b_size) + col] = block_c[row][col];
        }
      }
    }
  }

  std::vector<double> local_flat = FlattenMatrix(local_result);
  std::vector<int> recv_counts(size);
  std::vector<int> displs(size);

  int offset = 0;
  for (int proc = 0; proc < size; ++proc) {
    const int proc_block_rows = base_block_rows + (proc < extra_block_rows ? 1 : 0);
    recv_counts[proc] = proc_block_rows * b_size * n;
    displs[proc] = offset;
    offset += recv_counts[proc];
  }

  std::vector<double> global_flat(total_elems);
  MPI_Allgatherv(local_flat.data(), static_cast<int>(local_flat.size()), MPI_DOUBLE, global_flat.data(),
                 recv_counts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);

  GetOutput() = UnflattenMatrix(global_flat, n, n);
  return true;
}

bool TimurACannonMatrixMultiplicationALL::PostProcessingImpl() {
  return true;
}

}  // namespace timur_a_cannon
