#include "luchnikov_e_mult_of_dense_matrices_fox_algorithm_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <string>

#include "util/include/util.hpp"

namespace luchnikov_e_mult_of_dense_matrices_fox_algorithm_seq {

namespace {

bool LoadMatrixFromFile(const std::string &file_path, DenseMatrix &mat) {
  std::ifstream input(file_path);
  if (!input.is_open()) {
    return false;
  }

  input >> mat.rows >> mat.cols;
  if (input.fail() || mat.rows <= 0 || mat.cols <= 0) {
    return false;
  }

  std::size_t total = static_cast<std::size_t>(mat.rows) * mat.cols;
  mat.values.resize(total);
  for (std::size_t i = 0; i < total; ++i) {
    input >> mat.values[i];
    if (input.fail()) {
      return false;
    }
  }
  return !input.fail();
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

  // FIX: переменная 's' переименована в 'stage' для соответствия clang-tidy
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

  std::string task_id = "luchnikov_e_mult_of_dense_matrices_fox_algorithm_seq";
  path_a_ = ppc::util::GetAbsoluteTaskPath(task_id, "A_" + std::to_string(in) + ".txt");
  path_b_ = ppc::util::GetAbsoluteTaskPath(task_id, "B_" + std::to_string(in) + ".txt");
}

bool LuchnikovEMultOfDenseMatrixFoxAlgoritmSeq::ValidationImpl() {
  return GetInput() > 0 && GetOutput() == 0.0;
}

bool LuchnikovEMultOfDenseMatrixFoxAlgoritmSeq::PreProcessingImpl() {
  if (!LoadMatrixFromFile(path_a_, matrix_a_)) {
    return false;
  }
  if (!LoadMatrixFromFile(path_b_, matrix_b_)) {
    return false;
  }

  if (matrix_a_.rows != matrix_b_.rows || !matrix_a_.IsSquare()) {
    block_size_ = 0;
  } else {
    block_size_ = DetermineBlockSize(matrix_a_.rows);
  }
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
