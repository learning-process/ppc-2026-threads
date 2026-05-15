#include "yakimov_i_mult_of_dense_matrices_fox_algorithm_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <string>
#include <vector>

#include "util/include/util.hpp"
#include "yakimov_i_mult_of_dense_matrices_fox_algorithm_seq/common/include/common.hpp"

namespace yakimov_i_mult_of_dense_matrices_fox_algorithm_seq {

namespace {

bool ReadDimensions(std::ifstream &file, DenseMatrix &matrix) {
  bool success = true;

  file >> matrix.rows;
  file >> matrix.cols;

  if (matrix.rows <= 0 || matrix.cols <= 0) {
    success = false;
  }

  return success;
}

bool ReadMatrixData(std::ifstream &file, DenseMatrix &matrix) {
  bool success = true;

  auto total_elements = static_cast<std::size_t>(matrix.rows) * static_cast<std::size_t>(matrix.cols);
  matrix.data.resize(total_elements, 0.0);

  for (int i = 0; i < matrix.rows; ++i) {
    for (int j = 0; j < matrix.cols; ++j) {
      file >> matrix(i, j);
    }
  }

  return success;
}

bool ReadMatrixFromFileImpl(const std::string &filename, DenseMatrix &matrix) {
  std::ifstream file(filename);

  bool success = ReadDimensions(file, matrix);

  success = success && ReadMatrixData(file, matrix);

  file.close();

  return success;
}

void MultiplyBlock(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result, int row_start, int col_start,
                   int block_size, int a_row_offset, int b_col_offset) {
  for (int i = 0; i < block_size; ++i) {
    for (int j = 0; j < block_size; ++j) {
      double sum = 0.0;

      for (int k = 0; k < block_size; ++k) {
        double a_val = a(row_start + i, a_row_offset + k);
        double b_val = b(b_col_offset + k, col_start + j);

        sum += a_val * b_val;
      }

      result(row_start + i, col_start + j) += sum;
    }
  }
}

void FoxAlgorithmImpl(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result, int block_size) {
  if (block_size <= 0 || a.rows < block_size || a.cols < block_size || b.rows < block_size || b.cols < block_size) {
    return;
  }

  int num_blocks = a.rows / block_size;

  result.rows = a.rows;
  result.cols = b.cols;

  auto total_elements = static_cast<std::size_t>(result.rows) * static_cast<std::size_t>(result.cols);
  result.data.assign(total_elements, 0.0);

  for (int stage = 0; stage < num_blocks; ++stage) {
    for (int i = 0; i < num_blocks; ++i) {
      int broadcast_block = (i + stage) % num_blocks;

      for (int j = 0; j < num_blocks; ++j) {
        MultiplyBlock(a, b, result, i * block_size, j * block_size, block_size, broadcast_block * block_size,
                      j * block_size);
      }
    }
  }
}

}  // namespace

YakimovIMultOfDenseMatricesFoxAlgorithmSEQ::YakimovIMultOfDenseMatricesFoxAlgorithmSEQ(const InType &in) {
  this->SetTypeOfTask(YakimovIMultOfDenseMatricesFoxAlgorithmSEQ::GetStaticTypeOfTask());

  this->GetInput() = in;
  this->GetOutput() = 0.0;

  std::string task_name = "yakimov_i_mult_of_dense_matrices_fox_algorithm_seq";

  this->matrix_a_filename_ =
      ppc::util::GetAbsoluteTaskPath(task_name, "A_" + std::to_string(this->GetInput()) + ".txt");

  this->matrix_b_filename_ =
      ppc::util::GetAbsoluteTaskPath(task_name, "B_" + std::to_string(this->GetInput()) + ".txt");
}

bool YakimovIMultOfDenseMatricesFoxAlgorithmSEQ::ValidationImpl() {
  bool input_valid = (this->GetInput() > 0);

  bool output_valid = (this->GetOutput() == 0.0);

  return input_valid && output_valid;
}

bool YakimovIMultOfDenseMatricesFoxAlgorithmSEQ::PreProcessingImpl() {
  bool success = true;

  success = success && ReadMatrixFromFileImpl(this->matrix_a_filename_, this->matrix_a_);

  success = success && ReadMatrixFromFileImpl(this->matrix_b_filename_, this->matrix_b_);

  if (!success) {
    return false;
  }

  if (this->matrix_a_.cols != this->matrix_b_.rows) {
    return false;
  }

  int min_dimension = std::min(this->matrix_a_.rows, this->matrix_a_.cols);
  min_dimension = std::min(min_dimension, this->matrix_b_.rows);
  min_dimension = std::min(min_dimension, this->matrix_b_.cols);

  this->block_size_ = 1;
  while (this->block_size_ * 2 <= min_dimension) {
    this->block_size_ *= 2;
  }

  return this->block_size_ > 0;
}

bool YakimovIMultOfDenseMatricesFoxAlgorithmSEQ::RunImpl() {
  FoxAlgorithmImpl(this->matrix_a_, this->matrix_b_, this->result_matrix_, this->block_size_);

  return true;
}

bool YakimovIMultOfDenseMatricesFoxAlgorithmSEQ::PostProcessingImpl() {
  double sum = 0.0;

  for (double value : this->result_matrix_.data) {
    sum += value;
  }

  this->GetOutput() = sum;

  return true;
}

}  // namespace yakimov_i_mult_of_dense_matrices_fox_algorithm_seq
