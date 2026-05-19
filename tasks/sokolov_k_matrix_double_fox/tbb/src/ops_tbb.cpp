#include "sokolov_k_matrix_double_fox/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "oneapi/tbb/parallel_for.h"
#include "sokolov_k_matrix_double_fox/common/include/common.hpp"

namespace sokolov_k_matrix_double_fox {

namespace {

void DecomposeToBlocksTbb(const std::vector<double> &flat, std::vector<double> &blocks, int n, int bs, int q) {
  tbb::parallel_for(0, q, [&](int bi) {
    for (int bj = 0; bj < q; bj++) {
      int block_off = ((bi * q) + bj) * (bs * bs);
      for (int i = 0; i < bs; i++) {
        for (int j = 0; j < bs; j++) {
          blocks[block_off + (i * bs) + j] = flat[(((bi * bs) + i) * n) + ((bj * bs) + j)];
        }
      }
    }
  });
}

void AssembleFromBlocksTbb(const std::vector<double> &blocks, std::vector<double> &flat, int n, int bs, int q) {
  tbb::parallel_for(0, q, [&](int bi) {
    for (int bj = 0; bj < q; bj++) {
      int block_off = ((bi * q) + bj) * (bs * bs);
      for (int i = 0; i < bs; i++) {
        for (int j = 0; j < bs; j++) {
          flat[(((bi * bs) + i) * n) + ((bj * bs) + j)] = blocks[block_off + (i * bs) + j];
        }
      }
    }
  });
}

void MultiplyBlocksTbb(const std::vector<double> &a, int a_off, const std::vector<double> &b, int b_off,
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

void FoxStepTbb(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, int bs, int q,
                int step) {
  int bsq = bs * bs;
  tbb::parallel_for(0, q, [&](int i) {
    int k = (i + step) % q;
    for (int j = 0; j < q; j++) {
      MultiplyBlocksTbb(a, ((i * q) + k) * bsq, b, ((k * q) + j) * bsq, c, ((i * q) + j) * bsq, bs);
    }
  });
}

void FoxMultiplyTbb(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, int bs, int q) {
  for (int step = 0; step < q; step++) {
    FoxStepTbb(a, b, c, bs, q, step);
  }
}

int ChooseBlockSizeTbb(int n) {
  for (int div = static_cast<int>(std::sqrt(static_cast<double>(n))); div >= 1; div--) {
    if (n % div == 0) {
      return div;
    }
  }
  return 1;
}

}  // namespace

SokolovKMatrixDoubleFoxTBB::SokolovKMatrixDoubleFoxTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool SokolovKMatrixDoubleFoxTBB::ValidationImpl() {
  return (GetInput() > 0) && (GetOutput() == 0);
}

bool SokolovKMatrixDoubleFoxTBB::PreProcessingImpl() {
  GetOutput() = 0;
  n_ = GetInput();
  block_size_ = ChooseBlockSizeTbb(n_);
  q_ = n_ / block_size_;
  auto sz = static_cast<std::size_t>(n_) * n_;
  std::vector<double> a(sz, 1.5);
  std::vector<double> b(sz, 2.0);
  blocks_a_.resize(sz);
  blocks_b_.resize(sz);
  blocks_c_.assign(sz, 0.0);
  DecomposeToBlocksTbb(a, blocks_a_, n_, block_size_, q_);
  DecomposeToBlocksTbb(b, blocks_b_, n_, block_size_, q_);
  return true;
}

bool SokolovKMatrixDoubleFoxTBB::RunImpl() {
  std::ranges::fill(blocks_c_, 0.0);
  FoxMultiplyTbb(blocks_a_, blocks_b_, blocks_c_, block_size_, q_);
  return true;
}

bool SokolovKMatrixDoubleFoxTBB::PostProcessingImpl() {
  std::vector<double> result(static_cast<std::size_t>(n_) * n_);
  AssembleFromBlocksTbb(blocks_c_, result, n_, block_size_, q_);
  double expected = 3.0 * n_;
  bool ok = std::ranges::all_of(result, [expected](double v) { return std::abs(v - expected) <= 1e-9; });
  GetOutput() = ok ? GetInput() : -1;
  std::vector<double>().swap(blocks_a_);
  std::vector<double>().swap(blocks_b_);
  std::vector<double>().swap(blocks_c_);
  return true;
}

}  // namespace sokolov_k_matrix_double_fox
