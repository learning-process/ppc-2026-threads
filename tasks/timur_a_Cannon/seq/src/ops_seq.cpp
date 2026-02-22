#include "timur_a_Cannon/seq/include/ops_seq.hpp"

#include <cmath>

#include "timur_a_Cannon/common/include/common.hpp"

namespace timur_a_Cannon {

namespace {

void MultiplyAddBlock(const Matrix &A, const Matrix &B, Matrix &C, std::size_t i0, std::size_t j0, std::size_t k0,
                      std::size_t block_size) {
  const std::size_t n = A.n;
  for (std::size_t i = 0; i < block_size; ++i) {
    for (std::size_t k = 0; k < block_size; ++k) {
      double a_val = A(i0 + i, k0 + k);
      for (std::size_t j = 0; j < block_size; ++j) {
        C(i0 + i, j0 + j) += a_val * B(k0 + k, j0 + j);
      }
    }
  }
}

Matrix CannonMultiply(const Matrix &A, const Matrix &B) {
  std::size_t n = A.n;
  Matrix C(n, 0.0);

  if (n == 0) {
    return C;
  }

  std::size_t p = static_cast<std::size_t>(std::floor(std::sqrt(static_cast<double>(n))));
  while (p > 1 && (n % p != 0)) {
    --p;
  }

  if (p <= 1) {
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t k = 0; k < n; ++k) {
        double a_val = A(i, k);
        for (std::size_t j = 0; j < n; ++j) {
          C(i, j) += a_val * B(k, j);
        }
      }
    }
    return C;
  }

  const std::size_t block_size = n / p;

  for (std::size_t s = 0; s < p; ++s) {
    for (std::size_t bi = 0; bi < p; ++bi) {
      for (std::size_t bj = 0; bj < p; ++bj) {
        std::size_t a_block_col = (bj + bi + s) % p;
        std::size_t b_block_row = (bi + bj + s) % p;

        std::size_t i0 = bi * block_size;
        std::size_t j0 = bj * block_size;
        std::size_t k0 = a_block_col * block_size;

        MultiplyAddBlock(A, B, C, i0, j0, k0, block_size);
      }
    }
  }

  return C;
}

}  // namespace

TimurACannonSEQ::TimurACannonSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TimurACannonSEQ::ValidationImpl() {
  const auto &in = GetInput();
  if (in.A.n == 0 || in.B.n == 0) {
    return false;
  }
  if (in.A.n != in.B.n) {
    return false;
  }
  if (in.A.data.size() != in.A.n * in.A.n) {
    return false;
  }
  if (in.B.data.size() != in.B.n * in.B.n) {
    return false;
  }
  return true;
}

bool TimurACannonSEQ::PreProcessingImpl() {
  const auto &in = GetInput();
  GetOutput() = Matrix(in.A.n, 0.0);
  return true;
}

bool TimurACannonSEQ::RunImpl() {
  const auto &in = GetInput();
  GetOutput() = CannonMultiply(in.A, in.B);
  return true;
}

bool TimurACannonSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace timur_a_Cannon
