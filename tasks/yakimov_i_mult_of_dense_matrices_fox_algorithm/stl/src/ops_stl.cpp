#include "yakimov_i_mult_of_dense_matrices_fox_algorithm/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "util/include/util.hpp"
#include "yakimov_i_mult_of_dense_matrices_fox_algorithm/common/include/common.hpp"

namespace yakimov_i_mult_of_dense_matrices_fox_algorithm {

namespace {

bool ReadDimensions(std::ifstream &file, DenseMatrix &matrix) {
  file >> matrix.rows;
  file >> matrix.cols;
  return matrix.rows > 0 && matrix.cols > 0;
}

bool ReadMatrixData(std::ifstream &file, DenseMatrix &matrix) {
  auto total_elements = static_cast<std::size_t>(matrix.rows) * static_cast<std::size_t>(matrix.cols);
  matrix.data.resize(total_elements, 0.0);

  for (int i = 0; i < matrix.rows; ++i) {
    for (int j = 0; j < matrix.cols; ++j) {
      if (!(file >> matrix(i, j))) {
        return false;
      }
    }
  }
  return true;
}

bool ReadMatrixFromFileImpl(const std::string &filename, DenseMatrix &matrix) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    return false;
  }

  if (!ReadDimensions(file, matrix)) {
    return false;
  }
  if (!ReadMatrixData(file, matrix)) {
    return false;
  }

  file.close();
  return true;
}

void SimpleMultiply(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result) {
  result.rows = a.rows;
  result.cols = b.cols;
  result.data.assign(static_cast<std::size_t>(result.rows) * result.cols, 0.0);

  for (int i = 0; i < a.rows; ++i) {
    for (int j = 0; j < b.cols; ++j) {
      double sum = 0.0;
      for (int k = 0; k < a.cols; ++k) {
        sum += a(i, k) * b(k, j);
      }
      result(i, j) = sum;
    }
  }
}

void MultiplyBlock(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result, int row_start, int col_start,
                   int block_size, int a_row_offset, int b_col_offset) {
  for (int i = 0; i < block_size; ++i) {
    for (int j = 0; j < block_size; ++j) {
      double sum = 0.0;
      for (int k = 0; k < block_size; ++k) {
        sum += a(row_start + i, a_row_offset + k) * b(b_col_offset + k, col_start + j);
      }
      result(row_start + i, col_start + j) += sum;
    }
  }
}

void ProcessStage(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result, int stage, int num_blocks,
                  int block_size) {
  for (int i = 0; i < num_blocks; ++i) {
    int broadcast_block = (i + stage) % num_blocks;
    for (int j = 0; j < num_blocks; ++j) {
      MultiplyBlock(a, b, result, i * block_size, j * block_size, block_size, broadcast_block * block_size,
                    j * block_size);
    }
  }
}

unsigned int CalculateNumThreads(int num_blocks) {
  unsigned int hardware_threads = std::thread::hardware_concurrency();
  unsigned int num_threads = (hardware_threads == 0U) ? 2U : hardware_threads;
  unsigned int result = std::min(static_cast<unsigned int>(num_blocks), num_threads);
  return result;
}

void LaunchThreads(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result, int num_blocks, int block_size,
                   unsigned int num_threads) {
  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  int stages_per_thread = num_blocks / static_cast<int>(num_threads);
  int remaining_stages = num_blocks % static_cast<int>(num_threads);
  int stage_start = 0;

  for (unsigned int thread_index = 0U; thread_index < num_threads; ++thread_index) {
    int stages_for_this_thread = stages_per_thread;
    if (std::cmp_less(thread_index, remaining_stages)) {
      stages_for_this_thread = stages_for_this_thread + 1;
    }
    if (stages_for_this_thread == 0) {
      continue;
    }

    int stage_end = stage_start + stages_for_this_thread;

    threads.emplace_back([&a, &b, &result, stage_start, stage_end, num_blocks, block_size]() {
      for (int stage = stage_start; stage < stage_end; ++stage) {
        ProcessStage(a, b, result, stage, num_blocks, block_size);
      }
    });

    stage_start = stage_end;
  }

  for (std::thread &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

void FoxAlgorithmImpl(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result, int block_size) {
  if (block_size <= 0) {
    SimpleMultiply(a, b, result);
    return;
  }

  bool is_not_square = (a.rows != a.cols) || (b.rows != b.cols) || (a.rows != b.rows);
  if (is_not_square) {
    SimpleMultiply(a, b, result);
    return;
  }

  if (a.rows % block_size != 0) {
    SimpleMultiply(a, b, result);
    return;
  }

  int n = a.rows;
  int num_blocks = n / block_size;

  result.rows = n;
  result.cols = n;
  auto total_elements = static_cast<std::size_t>(n) * static_cast<std::size_t>(n);
  result.data.assign(total_elements, 0.0);

  unsigned int num_threads = CalculateNumThreads(num_blocks);
  LaunchThreads(a, b, result, num_blocks, block_size, num_threads);
}

}  // namespace

YakimovIMultOfDenseMatricesFoxAlgorithmSTL::YakimovIMultOfDenseMatricesFoxAlgorithmSTL(const InType &in) {
  this->SetTypeOfTask(YakimovIMultOfDenseMatricesFoxAlgorithmSTL::GetStaticTypeOfTask());
  this->GetInput() = in;
  this->GetOutput() = 0.0;

  std::string task_name = "yakimov_i_mult_of_dense_matrices_fox_algorithm";
  this->matrix_a_filename_ = ppc::util::GetAbsoluteTaskPath(task_name, "A_" + std::to_string(in) + ".txt");
  this->matrix_b_filename_ = ppc::util::GetAbsoluteTaskPath(task_name, "B_" + std::to_string(in) + ".txt");
}

bool YakimovIMultOfDenseMatricesFoxAlgorithmSTL::ValidationImpl() {
  return (this->GetInput() > 0) && (this->GetOutput() == 0.0);
}

bool YakimovIMultOfDenseMatricesFoxAlgorithmSTL::PreProcessingImpl() {
  if (!ReadMatrixFromFileImpl(this->matrix_a_filename_, this->matrix_a_)) {
    return false;
  }
  if (!ReadMatrixFromFileImpl(this->matrix_b_filename_, this->matrix_b_)) {
    return false;
  }
  if (this->matrix_a_.cols != this->matrix_b_.rows) {
    return false;
  }

  if (this->matrix_a_.rows != this->matrix_a_.cols || this->matrix_b_.rows != this->matrix_b_.cols ||
      this->matrix_a_.rows != this->matrix_b_.rows) {
    this->block_size_ = 0;
    return true;
  }

  int n = this->matrix_a_.rows;
  this->block_size_ = 64;
  while (this->block_size_ * 2 <= n && this->block_size_ < 256) {
    this->block_size_ *= 2;
  }

  this->block_size_ = std::min(this->block_size_, n);

  return this->block_size_ > 0;
}

bool YakimovIMultOfDenseMatricesFoxAlgorithmSTL::RunImpl() {
  FoxAlgorithmImpl(this->matrix_a_, this->matrix_b_, this->result_matrix_, this->block_size_);
  return true;
}

bool YakimovIMultOfDenseMatricesFoxAlgorithmSTL::PostProcessingImpl() {
  double sum = 0.0;
  for (double val : this->result_matrix_.data) {
    sum += val;
  }
  this->GetOutput() = sum;
  return true;
}

}  // namespace yakimov_i_mult_of_dense_matrices_fox_algorithm
