#include "yakimov_i_mult_of_dense_matrices_fox_algorithm_seq/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

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
      if (!(file >> matrix(i, j))) {
        success = false;
        break;
      }
    }
    if (!success) {
      break;
    }
  }

  return success;
}

bool ReadMatrixFromFileImpl(const std::string &filename, DenseMatrix &matrix) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    return false;
  }

  bool success = ReadDimensions(file, matrix);
  if (success) {
    success = ReadMatrixData(file, matrix);
  }

  file.close();
  return success;
}

void MultiplyBlockAtomic(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result, int row_start, int col_start,
                         int block_size, int a_row_offset, int b_col_offset) {
  for (int i = 0; i < block_size; ++i) {
    for (int j = 0; j < block_size; ++j) {
      double sum = 0.0;
      for (int k = 0; k < block_size; ++k) {
        double a_val = a(row_start + i, a_row_offset + k);
        double b_val = b(b_col_offset + k, col_start + j);
        sum += a_val * b_val;
      }
#pragma omp atomic
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

  const auto &a_local = a;
  const auto &b_local = b;
  auto &result_local = result;
  int num_blocks_local = num_blocks;
  int block_size_local = block_size;

#pragma omp parallel for collapse(2) default(none) \
    shared(a_local, b_local, result_local, num_blocks_local, block_size_local)
  for (int stage = 0; stage < num_blocks_local; ++stage) {
    for (int i = 0; i < num_blocks_local; ++i) {
      int broadcast_block = (i + stage) % num_blocks_local;
      for (int j = 0; j < num_blocks_local; ++j) {
        MultiplyBlockAtomic(a_local, b_local, result_local, i * block_size_local, j * block_size_local,
                            block_size_local, broadcast_block * block_size_local, j * block_size_local);
      }
    }
  }
}

}  // namespace

YakimovIMultOfDenseMatricesFoxAlgorithmOMP::YakimovIMultOfDenseMatricesFoxAlgorithmOMP(const InType &in) {
  this->SetTypeOfTask(YakimovIMultOfDenseMatricesFoxAlgorithmOMP::GetStaticTypeOfTask());

  this->GetInput() = in;
  this->GetOutput() = 0.0;

  std::string task_name = "yakimov_i_mult_of_dense_matrices_fox_algorithm_seq";

  this->matrix_a_filename_ =
      ppc::util::GetAbsoluteTaskPath(task_name, "A_" + std::to_string(this->GetInput()) + ".txt");

  this->matrix_b_filename_ =
      ppc::util::GetAbsoluteTaskPath(task_name, "B_" + std::to_string(this->GetInput()) + ".txt");

#ifdef _OPENMP
  // Инициализация OpenMP
  omp_set_num_threads(omp_get_max_threads());
#endif
}

YakimovIMultOfDenseMatricesFoxAlgorithmOMP::~YakimovIMultOfDenseMatricesFoxAlgorithmOMP() {
#ifdef _OPENMP
// Синхронизация всех потоков перед завершением
#  pragma omp parallel
  {
#  pragma omp single
    {
      // Пустая директива для синхронизации
    }
  }
#endif
}

bool YakimovIMultOfDenseMatricesFoxAlgorithmOMP::ValidationImpl() {
  bool input_valid = (this->GetInput() > 0);
  bool output_valid = (this->GetOutput() == 0.0);
  return input_valid && output_valid;
}

bool YakimovIMultOfDenseMatricesFoxAlgorithmOMP::PreProcessingImpl() {
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

  this->block_size_ = 64;
  while (this->block_size_ * 2 <= min_dimension && this->block_size_ < 256) {
    this->block_size_ *= 2;
  }

  return this->block_size_ > 0;
}

bool YakimovIMultOfDenseMatricesFoxAlgorithmOMP::RunImpl() {
  FoxAlgorithmImpl(this->matrix_a_, this->matrix_b_, this->result_matrix_, this->block_size_);
  return true;
}

bool YakimovIMultOfDenseMatricesFoxAlgorithmOMP::PostProcessingImpl() {
  double sum = 0.0;

  for (double val : this->result_matrix_.data) {
    sum += val;
  }

  this->GetOutput() = sum;

  this->matrix_a_.data.clear();
  this->matrix_a_.data.shrink_to_fit();
  this->matrix_b_.data.clear();
  this->matrix_b_.data.shrink_to_fit();
  this->result_matrix_.data.clear();
  this->result_matrix_.data.shrink_to_fit();

  return true;
}

}  // namespace yakimov_i_mult_of_dense_matrices_fox_algorithm_seq
