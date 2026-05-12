#include "yakimov_i_mult_of_dense_matrices_fox_algorithm/all/include/ops_all.hpp"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <tbb/parallel_for.h>
#include <tbb/task_group.h>

#include "util/include/util.hpp"
#include "yakimov_i_mult_of_dense_matrices_fox_algorithm/common/include/common.hpp"

namespace yakimov_i_mult_of_dense_matrices_fox_algorithm {

namespace {

constexpr int kSmallMatrixThreshold = 64;
constexpr int kMediumMatrixThreshold = 256;
constexpr int kBlockSizeSmall = 32;
constexpr int kBlockSizeMedium = 64;
constexpr int kBlockSizeLarge = 128;
constexpr int kUnrollFactor = 4;

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

void SimpleMultiplySeq(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result) {
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

void SimpleMultiplyOMP(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result) {
  result.rows = a.rows;
  result.cols = b.cols;
  result.data.assign(static_cast<std::size_t>(result.rows) * result.cols, 0.0);

#ifdef _OPENMP
#  pragma omp parallel for schedule(static) default(none) shared(a, b, result)
#endif
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

void SimpleMultiplyTBB(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result) {
  result.rows = a.rows;
  result.cols = b.cols;
  result.data.assign(static_cast<std::size_t>(result.rows) * result.cols, 0.0);

  tbb::parallel_for(0, a.rows, [&](int i) {
    for (int j = 0; j < b.cols; ++j) {
      double sum = 0.0;
      for (int k = 0; k < a.cols; ++k) {
        sum += a(i, k) * b(k, j);
      }
      result(i, j) = sum;
    }
  });
}

void MultiplyBlockUnrolled(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result, int row_start,
                           int col_start, int block_size, int a_row_offset, int b_col_offset) {
  int block_size_aligned = (block_size / kUnrollFactor) * kUnrollFactor;

  for (int i = 0; i < block_size; ++i) {
    for (int j = 0; j < block_size; ++j) {
      double sum = 0.0;
      int k = 0;

      for (; k < block_size_aligned; k += kUnrollFactor) {
        sum += a(row_start + i, a_row_offset + k) * b(b_col_offset + k, col_start + j);
        sum += a(row_start + i, a_row_offset + k + 1) * b(b_col_offset + k + 1, col_start + j);
        sum += a(row_start + i, a_row_offset + k + 2) * b(b_col_offset + k + 2, col_start + j);
        sum += a(row_start + i, a_row_offset + k + 3) * b(b_col_offset + k + 3, col_start + j);
      }

      for (; k < block_size; ++k) {
        sum += a(row_start + i, a_row_offset + k) * b(b_col_offset + k, col_start + j);
      }

      result(row_start + i, col_start + j) += sum;
    }
  }
}

void ProcessStageOMP(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result, int stage, int num_blocks,
                     int block_size) {
#ifdef _OPENMP
#  pragma omp parallel for schedule(dynamic) default(none) shared(a, b, result, stage, num_blocks, block_size)
#endif
  for (int idx = 0; idx < num_blocks * num_blocks; ++idx) {
    int i = idx / num_blocks;
    int j = idx % num_blocks;
    int broadcast_block = (i + stage) % num_blocks;
    MultiplyBlockUnrolled(a, b, result, i * block_size, j * block_size, block_size, broadcast_block * block_size,
                          j * block_size);
  }
}

void ProcessStageTBB(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result, int stage, int num_blocks,
                     int block_size) {
  for (int i = 0; i < num_blocks; ++i) {
    int broadcast_block = (i + stage) % num_blocks;
    for (int j = 0; j < num_blocks; ++j) {
      MultiplyBlockUnrolled(a, b, result, i * block_size, j * block_size, block_size, broadcast_block * block_size,
                            j * block_size);
    }
  }
}

void HandleNonSquareMatrices(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result) {
  int n = a.rows;

  if (n <= kSmallMatrixThreshold) {
    SimpleMultiplySeq(a, b, result);
  } else if (n <= kMediumMatrixThreshold) {
    SimpleMultiplyOMP(a, b, result);
  } else {
    SimpleMultiplyTBB(a, b, result);
  }
}

void HandleSmallNumBlocks(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result, int num_blocks,
                          int block_size) {
  for (int stage = 0; stage < num_blocks; ++stage) {
    ProcessStageOMP(a, b, result, stage, num_blocks, block_size);
  }
}

void HandleMediumNumBlocks(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result, int num_blocks,
                           int block_size) {
  tbb::parallel_for(0, num_blocks, [&](int stage) { ProcessStageOMP(a, b, result, stage, num_blocks, block_size); });
}

void HandleLargeNumBlocks(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result, int num_blocks,
                          int block_size) {
  unsigned int hardware_threads = std::thread::hardware_concurrency();
  int num_tasks = static_cast<int>(hardware_threads);
  num_tasks = std::max(2, std::min(num_tasks, num_blocks));
  int stages_per_task = num_blocks / num_tasks;
  int remaining_stages = num_blocks % num_tasks;

  tbb::task_group task_group;

  int stage_start = 0;
  for (int task_index = 0; task_index < num_tasks; ++task_index) {
    int stages_for_this_task = stages_per_task;
    if (task_index < remaining_stages) {
      stages_for_this_task = stages_for_this_task + 1;
    }
    if (stages_for_this_task == 0) {
      continue;
    }

    int stage_end = stage_start + stages_for_this_task;

    task_group.run([&a, &b, &result, stage_start, stage_end, num_blocks, block_size]() {
      for (int stage = stage_start; stage < stage_end; ++stage) {
        ProcessStageTBB(a, b, result, stage, num_blocks, block_size);
      }
    });

    stage_start = stage_end;
  }

  task_group.wait();
}

void FoxAlgorithmAdaptive(const DenseMatrix &a, const DenseMatrix &b, DenseMatrix &result, int block_size) {
  bool is_not_square = (a.rows != a.cols) || (b.rows != b.cols) || (a.rows != b.rows);
  bool is_not_divisible = (a.rows % block_size != 0);

  if (is_not_square || is_not_divisible) {
    HandleNonSquareMatrices(a, b, result);
    return;
  }

  int n = a.rows;
  int num_blocks = n / block_size;
  result.rows = n;
  result.cols = n;
  result.data.assign(static_cast<std::size_t>(n) * n, 0.0);

  if (num_blocks <= 8) {
    HandleSmallNumBlocks(a, b, result, num_blocks, block_size);
  } else if (num_blocks <= 32) {
    HandleMediumNumBlocks(a, b, result, num_blocks, block_size);
  } else {
    HandleLargeNumBlocks(a, b, result, num_blocks, block_size);
  }
}

}  // namespace

YakimovIMultOfDenseMatricesFoxAlgorithmAll::YakimovIMultOfDenseMatricesFoxAlgorithmAll(const InType &in) {
  this->SetTypeOfTask(YakimovIMultOfDenseMatricesFoxAlgorithmAll::GetStaticTypeOfTask());
  this->GetInput() = in;
  this->GetOutput() = 0.0;

  std::string task_name = "yakimov_i_mult_of_dense_matrices_fox_algorithm";
  this->matrix_a_filename_ = ppc::util::GetAbsoluteTaskPath(task_name, "A_" + std::to_string(in) + ".txt");
  this->matrix_b_filename_ = ppc::util::GetAbsoluteTaskPath(task_name, "B_" + std::to_string(in) + ".txt");
}

bool YakimovIMultOfDenseMatricesFoxAlgorithmAll::ValidationImpl() {
  return (this->GetInput() > 0) && (this->GetOutput() == 0.0);
}

bool YakimovIMultOfDenseMatricesFoxAlgorithmAll::PreProcessingImpl() {
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

  if (n <= kSmallMatrixThreshold) {
    this->block_size_ = kBlockSizeSmall;
  } else if (n <= kMediumMatrixThreshold) {
    this->block_size_ = kBlockSizeMedium;
  } else {
    this->block_size_ = kBlockSizeLarge;
  }

  while (this->block_size_ * 2 <= n && this->block_size_ < kBlockSizeLarge) {
    this->block_size_ *= 2;
  }

  this->block_size_ = std::min(this->block_size_, n);

  return this->block_size_ > 0;
}

bool YakimovIMultOfDenseMatricesFoxAlgorithmAll::RunImpl() {
  FoxAlgorithmAdaptive(this->matrix_a_, this->matrix_b_, this->result_matrix_, this->block_size_);
  return true;
}

bool YakimovIMultOfDenseMatricesFoxAlgorithmAll::PostProcessingImpl() {
  double sum = 0.0;
  for (double value : this->result_matrix_.data) {
    sum += value;
  }
  this->GetOutput() = sum;
  return true;
}

}  // namespace yakimov_i_mult_of_dense_matrices_fox_algorithm
