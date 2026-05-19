#include "luchnikov_e_mult_of_dense_matrices_fox_algorithm_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>

namespace luchnikov_e_mult_of_dense_matrices_fox_algorithm_seq {

namespace {

void FillMatrixWithOnes(DenseMatrix &mat, int n) {
  mat.rows = n;
  mat.cols = n;
  mat.values.assign(static_cast<std::size_t>(n) * n, 1.0);
}

void FillMatrixWithZeros(DenseMatrix &mat, int r, int c) {
  mat.rows = r;
  mat.cols = c;
  mat.values.assign(static_cast<std::size_t>(r) * c, 0.0);
}

void NaiveMultiply(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &res) {
  FillMatrixWithZeros(res, a.rows, b.cols);
  for (int i = 0; i < a.rows; ++i) {
    for (int k = 0; k < a.cols; ++k) {
      double factor = a.At(i, k);
      if (factor == 0.0) {
        continue;
      }
      for (int j = 0; j < b.cols; ++j) {
        res.At(i, j) += factor * b.At(k, j);
      }
    }
  }
}

void MultiplyBlocks(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &res, int row_off, int col_off, int blk,
                    int a_col_shift, int b_row_shift) {
  for (int i = 0; i < blk; ++i) {
    for (int j = 0; j < blk; ++j) {
      double acc = 0.0;
      for (int k = 0; k < blk; ++k) {
        acc += a.At(row_off + i, a_col_shift + k) * b.At(b_row_shift + k, col_off + j);
      }
      res.At(row_off + i, col_off + j) += acc;
    }
  }
}

int DetermineBlockSize(int n) {
  int blk = n / 4;
  if (blk <= 0) {
    blk = 1;
  }
  return std::min(blk, 128);
}

void ExecuteFoxAlgorithm(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &res, int blk) {
  if (!a.IsSquare() || !b.IsSquare() || a.rows != b.rows || blk <= 0 || a.rows % blk != 0) {
    NaiveMultiply(a, b, res);
    return;
  }

  int n = a.rows;
  int stages = n / blk;
  FillMatrixWithZeros(res, n, n);

  for (int stage = 0; stage < stages; ++stage) {
    for (int i = 0; i < stages; ++i) {
      int broadcast_idx = (i + stage) % stages;
      for (int j = 0; j < stages; ++j) {
        MultiplyBlocks(a, b, res, i * blk, j * blk, blk, broadcast_idx * blk, broadcast_idx * blk);
      }
    }
  }
}

}  // namespace

LuchnikovEMultOfDenseMatrixFoxAlgoritmSeq::LuchnikovEMultOfDenseMatrixFoxAlgoritmSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool LuchnikovEMultOfDenseMatrixFoxAlgoritmSeq::ValidationImpl() {
  return GetInput() > 0 && GetOutput() == 0.0;
}

bool LuchnikovEMultOfDenseMatrixFoxAlgoritmSeq::PreProcessingImpl() {
  int n = GetInput();
  FillMatrixWithOnes(matrix_a_, n);
  FillMatrixWithOnes(matrix_b_, n);

  block_size_ = DetermineBlockSize(n);
  return true;
}

bool LuchnikovEMultOfDenseMatrixFoxAlgoritmSeq::RunImpl() {
  ExecuteFoxAlgorithm(matrix_a_, matrix_b_, result_matrix_, block_size_);
  return true;
}

bool LuchnikovEMultOfDenseMatrixFoxAlgoritmSeq::PostProcessingImpl() {
  double total = 0.0;
  for (double val : result_matrix_.values) {
    total += val;
  }
  GetOutput() = total;
  return std::isfinite(total);
}

}  // namespace luchnikov_e_mult_of_dense_matrices_fox_algorithm_seq
