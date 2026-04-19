#include "sokolov_k_matrix_double_fox/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "sokolov_k_matrix_double_fox/common/include/common.hpp"

namespace sokolov_k_matrix_double_fox {

namespace {

void DecomposeToBlocksAll(const std::vector<double> &flat, std::vector<double> &blocks, int n, int bs, int q) {
  for (int bi = 0; bi < q; bi++) {
    for (int bj = 0; bj < q; bj++) {
      int block_off = ((bi * q) + bj) * (bs * bs);
      for (int i = 0; i < bs; i++) {
        for (int j = 0; j < bs; j++) {
          blocks[block_off + (i * bs) + j] = flat[(((bi * bs) + i) * n) + ((bj * bs) + j)];
        }
      }
    }
  }
}

void AssembleFromBlocksAll(const std::vector<double> &blocks, std::vector<double> &flat, int n, int bs, int q) {
  for (int bi = 0; bi < q; bi++) {
    for (int bj = 0; bj < q; bj++) {
      int block_off = ((bi * q) + bj) * (bs * bs);
      for (int i = 0; i < bs; i++) {
        for (int j = 0; j < bs; j++) {
          flat[(((bi * bs) + i) * n) + ((bj * bs) + j)] = blocks[block_off + (i * bs) + j];
        }
      }
    }
  }
}

void MultiplyBlocksAll(const double *a, const double *b, double *c, int bs) {
  for (int i = 0; i < bs; i++) {
    for (int k = 0; k < bs; k++) {
      double val = a[(i * bs) + k];
      for (int j = 0; j < bs; j++) {
        c[(i * bs) + j] += val * b[(k * bs) + j];
      }
    }
  }
}

void FoxStepMpiOmp(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, int bs, int q,
                   int step, int row_begin, int row_end) {
  int bsq = bs * bs;
#pragma omp parallel for default(none) shared(a, b, c, bs, q, bsq, step, row_begin, row_end) schedule(static)
  for (int i = row_begin; i < row_end; i++) {
    int k = (i + step) % q;
    for (int j = 0; j < q; j++) {
      int a_off = ((i * q) + k) * bsq;
      int b_off = ((k * q) + j) * bsq;
      int c_off = ((i * q) + j) * bsq;
      MultiplyBlocksAll(a.data() + a_off, b.data() + b_off, c.data() + c_off, bs);
    }
  }
}

int ChooseBlockSizeAll(int n) {
  for (int div = static_cast<int>(std::sqrt(static_cast<double>(n))); div >= 1; div--) {
    if (n % div == 0) {
      return div;
    }
  }
  return 1;
}

void ComputeRowRange(int rank, int num_procs, int rows_per, int leftover, int &row_start, int &row_count) {
  if (rank < num_procs) {
    row_start = (rank * rows_per) + std::min(rank, leftover);
    row_count = rows_per + (rank < leftover ? 1 : 0);
  } else {
    row_start = 0;
    row_count = 0;
  }
}

void GatherResults(std::vector<double> &blocks_c, int rank, int num_procs, int rows_per, int leftover, int q, int bsq) {
  if (rank == 0) {
    for (int pr = 1; pr < num_procs; pr++) {
      int pr_start = 0;
      int pr_count = 0;
      ComputeRowRange(pr, num_procs, rows_per, leftover, pr_start, pr_count);
      if (pr_count > 0) {
        int offset = pr_start * q * bsq;
        int count = pr_count * q * bsq;
        MPI_Recv(blocks_c.data() + offset, count, MPI_DOUBLE, pr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
  } else if (rank < num_procs) {
    int my_start = 0;
    int my_count = 0;
    ComputeRowRange(rank, num_procs, rows_per, leftover, my_start, my_count);
    if (my_count > 0) {
      int offset = my_start * q * bsq;
      int count = my_count * q * bsq;
      MPI_Send(blocks_c.data() + offset, count, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
  }
}

}  // namespace

SokolovKMatrixDoubleFoxALL::SokolovKMatrixDoubleFoxALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool SokolovKMatrixDoubleFoxALL::ValidationImpl() {
  return (GetInput() > 0) && (GetOutput() == 0);
}

bool SokolovKMatrixDoubleFoxALL::PreProcessingImpl() {
  GetOutput() = 0;
  n_ = GetInput();
  block_size_ = ChooseBlockSizeAll(n_);
  q_ = n_ / block_size_;
  auto sz = static_cast<std::size_t>(n_) * n_;
  std::vector<double> a(sz, 1.5);
  std::vector<double> b(sz, 2.0);
  blocks_a_.resize(sz);
  blocks_b_.resize(sz);
  blocks_c_.assign(sz, 0.0);
  DecomposeToBlocksAll(a, blocks_a_, n_, block_size_, q_);
  DecomposeToBlocksAll(b, blocks_b_, n_, block_size_, q_);
  return true;
}

bool SokolovKMatrixDoubleFoxALL::RunImpl() {
  std::ranges::fill(blocks_c_, 0.0);

  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int total = static_cast<int>(blocks_a_.size());
  MPI_Bcast(&n_, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&block_size_, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&q_, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&total, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    blocks_a_.resize(total);
    blocks_b_.resize(total);
    blocks_c_.assign(total, 0.0);
  }

  MPI_Bcast(blocks_a_.data(), total, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(blocks_b_.data(), total, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  int num_procs = std::min(world_size, q_);
  int rows_per = q_ / std::max(num_procs, 1);
  int leftover = q_ % std::max(num_procs, 1);

  int my_row_start = 0;
  int my_row_count = 0;
  ComputeRowRange(rank, num_procs, rows_per, leftover, my_row_start, my_row_count);

  for (int step = 0; step < q_; step++) {
    FoxStepMpiOmp(blocks_a_, blocks_b_, blocks_c_, block_size_, q_, step, my_row_start, my_row_start + my_row_count);
  }

  int bsq = block_size_ * block_size_;
  GatherResults(blocks_c_, rank, num_procs, rows_per, leftover, q_, bsq);

  MPI_Bcast(blocks_c_.data(), total, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

bool SokolovKMatrixDoubleFoxALL::PostProcessingImpl() {
  std::vector<double> result(static_cast<std::size_t>(n_) * n_);
  AssembleFromBlocksAll(blocks_c_, result, n_, block_size_, q_);
  double expected = 3.0 * n_;
  bool ok = std::ranges::all_of(result, [expected](double v) { return std::abs(v - expected) <= 1e-9; });
  GetOutput() = ok ? GetInput() : -1;
  std::vector<double>().swap(blocks_a_);
  std::vector<double>().swap(blocks_b_);
  std::vector<double>().swap(blocks_c_);
  return true;
}

}  // namespace sokolov_k_matrix_double_fox
