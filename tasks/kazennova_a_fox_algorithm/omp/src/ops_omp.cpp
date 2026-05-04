#include "kazennova_a_fox_algorithm/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "kazennova_a_fox_algorithm/common/include/common.hpp"

namespace kazennova_a_fox_algorithm {

namespace {

int ChooseBlockSize(int n) {
  int best = 1;
  int sqrt_n = static_cast<int>(std::sqrt(static_cast<double>(n)));

  for (int bs = sqrt_n; bs >= 1; --bs) {
    if (n % bs == 0) {
      best = bs;
      break;
    }
  }
  return best;
}

}  // namespace

KazennovaATestTaskOMP::KazennovaATestTaskOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KazennovaATestTaskOMP::ValidationImpl() {
  const auto &in = GetInput();

  if (in.A.data.empty() || in.B.data.empty()) {
    return false;
  }

  if (in.A.rows <= 0 || in.A.cols <= 0 || in.B.rows <= 0 || in.B.cols <= 0) {
    return false;
  }

  if (in.A.cols != in.B.rows) {
    return false;
  }

  if (in.A.rows != in.A.cols || in.B.rows != in.B.cols || in.A.rows != in.B.rows) {
    return false;
  }

  return true;
}

void KazennovaATestTaskOMP::DecomposeMatrix(const std::vector<double> &src, std::vector<double> &dst, int n, int bs,
                                            int q) {
  size_t block_elements = static_cast<size_t>(bs) * bs;

  for (int bi = 0; bi < q; ++bi) {
    for (int bj = 0; bj < q; ++bj) {
      size_t block_offset = ((static_cast<size_t>(bi) * q) + bj) * block_elements;

      for (int i = 0; i < bs; ++i) {
        for (int j = 0; j < bs; ++j) {
          size_t src_idx = (((static_cast<size_t>(bi) * bs) + i) * n) + ((static_cast<size_t>(bj) * bs) + j);
          size_t dst_idx = block_offset + (static_cast<size_t>(i) * bs) + j;
          dst[dst_idx] = src[src_idx];
        }
      }
    }
  }
}

void KazennovaATestTaskOMP::AssembleMatrix(const std::vector<double> &src, std::vector<double> &dst, int n, int bs,
                                           int q) {
  size_t block_elements = static_cast<size_t>(bs) * bs;

  for (int bi = 0; bi < q; ++bi) {
    for (int bj = 0; bj < q; ++bj) {
      size_t block_offset = ((static_cast<size_t>(bi) * q) + bj) * block_elements;

      for (int i = 0; i < bs; ++i) {
        for (int j = 0; j < bs; ++j) {
          size_t dst_idx = (((static_cast<size_t>(bi) * bs) + i) * n) + ((static_cast<size_t>(bj) * bs) + j);
          size_t src_idx = block_offset + (static_cast<size_t>(i) * bs) + j;
          dst[dst_idx] = src[src_idx];
        }
      }
    }
  }
}

bool KazennovaATestTaskOMP::PreProcessingImpl() {
  const auto &in = GetInput();

  matrix_size_ = in.A.rows;
  GetOutput().rows = matrix_size_;
  GetOutput().cols = matrix_size_;
  GetOutput().data.assign(static_cast<size_t>(matrix_size_) * matrix_size_, 0.0);

  block_size_ = ChooseBlockSize(matrix_size_);
  block_count_ = matrix_size_ / block_size_;

  size_t total_blocks = static_cast<size_t>(block_count_) * block_count_;
  size_t block_elements = static_cast<size_t>(block_size_) * block_size_;
  a_blocks_.assign(total_blocks * block_elements, 0.0);
  b_blocks_.assign(total_blocks * block_elements, 0.0);
  c_blocks_.assign(total_blocks * block_elements, 0.0);

  DecomposeMatrix(in.A.data, a_blocks_, matrix_size_, block_size_, block_count_);
  DecomposeMatrix(in.B.data, b_blocks_, matrix_size_, block_size_, block_count_);

  return true;
}

void KazennovaATestTaskOMP::MultiplyBlock(size_t a_idx, size_t b_idx, size_t c_idx, int bs) {
  for (int ii = 0; ii < bs; ++ii) {
    size_t row_offset = a_idx + (static_cast<size_t>(ii) * bs);
    size_t c_row_offset = c_idx + (static_cast<size_t>(ii) * bs);
    for (int kk = 0; kk < bs; ++kk) {
      double a_val = a_blocks_[row_offset + kk];
      size_t b_row_offset = b_idx + (static_cast<size_t>(kk) * bs);
      for (int jj = 0; jj < bs; ++jj) {
        c_blocks_[c_row_offset + jj] += a_val * b_blocks_[b_row_offset + jj];
      }
    }
  }
}

bool KazennovaATestTaskOMP::RunImpl() {
  size_t block_elements = static_cast<size_t>(block_size_) * block_size_;

  for (auto &c_block : c_blocks_) {
    c_block = 0.0;
  }

  for (int step = 0; step < block_count_; ++step) {
#pragma omp parallel for default(none) shared(step, block_elements)
    for (int i = 0; i < block_count_; ++i) {
      for (int j = 0; j < block_count_; ++j) {
        int k = (i + step) % block_count_;

        size_t a_idx = ((static_cast<size_t>(i) * block_count_) + k) * block_elements;
        size_t b_idx = ((static_cast<size_t>(k) * block_count_) + j) * block_elements;
        size_t c_idx = ((static_cast<size_t>(i) * block_count_) + j) * block_elements;

        MultiplyBlock(a_idx, b_idx, c_idx, block_size_);
      }
    }
  }

  return true;
}

bool KazennovaATestTaskOMP::PostProcessingImpl() {
  AssembleMatrix(c_blocks_, GetOutput().data, matrix_size_, block_size_, block_count_);
  return true;
}

}  // namespace kazennova_a_fox_algorithm
