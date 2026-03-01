#include "sokolov_k_matrix_double_fox_seq/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

#include "sokolov_k_matrix_double_fox_seq/common/include/common.hpp"

namespace sokolov_k_matrix_double_fox_seq {

namespace {

void DecomposeToBlocks(const std::vector<double> &flat, std::vector<double> &blocks, int n, int bs, int q) {
  for (int bi = 0; bi < q; bi++) {
    for (int bj = 0; bj < q; bj++) {
      int block_off = (bi * q + bj) * bs * bs;
      for (int i = 0; i < bs; i++) {
        for (int j = 0; j < bs; j++) {
          blocks[block_off + (i * bs) + j] = flat[((bi * bs) + i) * n + (bj * bs + j)];
        }
      }
    }
  }
}

void AssembleFromBlocks(const std::vector<double> &blocks, std::vector<double> &flat, int n, int bs, int q) {
  for (int bi = 0; bi < q; bi++) {
    for (int bj = 0; bj < q; bj++) {
      int block_off = (bi * q + bj) * bs * bs;
      for (int i = 0; i < bs; i++) {
        for (int j = 0; j < bs; j++) {
          flat[((bi * bs) + i) * n + (bj * bs + j)] = blocks[block_off + (i * bs) + j];
        }
      }
    }
  }
}

void MultiplyBlocks(const std::vector<double> &a, int a_off, const std::vector<double> &b, int b_off,
                    std::vector<double> &c, int c_off, int bs) {
  for (int i = 0; i < bs; i++) {
    for (int k = 0; k < bs; k++) {
      double val = a[a_off + (i * bs) + k];
      for (int j = 0; j < bs; j++) {
        c[c_off + (i * bs) + j] += val * b[b_off + (k * bs) + j];
      }
    }
  }
}

void FoxStep(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, int bs, int q,
             int step) {
  int bsq = bs * bs;
  for (int i = 0; i < q; i++) {
    int k = (i + step) % q;
    for (int j = 0; j < q; j++) {
      MultiplyBlocks(a, (i * q + k) * bsq, b, (k * q + j) * bsq, c, (i * q + j) * bsq, bs);
    }
  }
}

void FoxMultiply(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, int bs, int q) {
  for (int step = 0; step < q; step++) {
    FoxStep(a, b, c, bs, q, step);
  }
}

}  // namespace

SokolovKMatrixDoubleFoxSEQ::SokolovKMatrixDoubleFoxSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SokolovKMatrixDoubleFoxSEQ::ValidationImpl() {
  const auto &in = GetInput();
  int n = std::get<0>(in);
  int b = std::get<1>(in);
  const auto &a = std::get<2>(in);
  const auto &mat_b = std::get<3>(in);
  if (n <= 0 || b <= 0 || (n % b != 0)) {
    return false;
  }
  auto expected = static_cast<std::size_t>(n) * n;
  if (a.size() != expected || mat_b.size() != expected) {
    return false;
  }
  for (std::size_t i = 0; i < a.size(); i++) {
    if (!std::isfinite(a[i]) || !std::isfinite(mat_b[i])) {
      return false;
    }
  }
  return true;
}

bool SokolovKMatrixDoubleFoxSEQ::PreProcessingImpl() {
  const auto &in = GetInput();
  n_ = std::get<0>(in);
  block_size_ = std::get<1>(in);
  q_ = n_ / block_size_;
  auto sz = static_cast<std::size_t>(n_) * n_;
  blocks_a_.resize(sz);
  blocks_b_.resize(sz);
  blocks_c_.assign(sz, 0.0);
  DecomposeToBlocks(std::get<2>(in), blocks_a_, n_, block_size_, q_);
  DecomposeToBlocks(std::get<3>(in), blocks_b_, n_, block_size_, q_);
  return true;
}

bool SokolovKMatrixDoubleFoxSEQ::RunImpl() {
  FoxMultiply(blocks_a_, blocks_b_, blocks_c_, block_size_, q_);
  return true;
}

bool SokolovKMatrixDoubleFoxSEQ::PostProcessingImpl() {
  GetOutput().resize(static_cast<std::size_t>(n_) * n_);
  AssembleFromBlocks(blocks_c_, GetOutput(), n_, block_size_, q_);
  std::vector<double>().swap(blocks_a_);
  std::vector<double>().swap(blocks_b_);
  std::vector<double>().swap(blocks_c_);
  return true;
}

}  // namespace sokolov_k_matrix_double_fox_seq
