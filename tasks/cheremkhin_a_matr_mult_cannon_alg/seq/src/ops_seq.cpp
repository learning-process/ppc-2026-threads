#include "cheremkhin_a_matr_mult_cannon_alg/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "cheremkhin_a_matr_mult_cannon_alg/common/include/common.hpp"

namespace cheremkhin_a_matr_mult_cannon_alg {

namespace {

inline std::size_t Idx(std::size_t n, std::size_t r, std::size_t c) {
  return (r * n) + c;
}

std::size_t ChooseQ(std::size_t n) {
  if (n <= 1) {
    return 1;
  }
  const auto root = static_cast<std::size_t>(std::sqrt(static_cast<double>(n)));
  return (root == 0) ? 1 : root;
}

std::size_t CeilDiv(std::size_t a, std::size_t b) {
  return (a + b - 1) / b;
}

void MulAddBlock(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, std::size_t n,
                 std::size_t bs, std::size_t bi, std::size_t bk, std::size_t bj) {
  const std::size_t i0 = bi * bs;
  const std::size_t k0 = bk * bs;
  const std::size_t j0 = bj * bs;

  for (std::size_t ii = 0; ii < bs; ++ii) {
    const std::size_t i = i0 + ii;
    const std::size_t a_row = i * n;
    const std::size_t c_row = i * n;
    for (std::size_t kk = 0; kk < bs; ++kk) {
      const std::size_t k = k0 + kk;
      const double aik = a[a_row + k];
      const std::size_t b_row = k * n;
      for (std::size_t jj = 0; jj < bs; ++jj) {
        const std::size_t j = j0 + jj;
        c[c_row + j] += aik * b[b_row + j];
      }
    }
  }
}

}  // namespace

CheremkhinAMatrMultCannonAlgSEQ::CheremkhinAMatrMultCannonAlgSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool CheremkhinAMatrMultCannonAlgSEQ::ValidationImpl() {
  const std::size_t n = std::get<0>(GetInput());
  const auto &a = std::get<1>(GetInput());
  const auto &b = std::get<2>(GetInput());
  return n > 0 && a.size() == n * n && b.size() == n * n;
}

bool CheremkhinAMatrMultCannonAlgSEQ::PreProcessingImpl() {
  GetOutput() = {};
  return true;
}

bool CheremkhinAMatrMultCannonAlgSEQ::RunImpl() {
  const std::size_t n = std::get<0>(GetInput());
  const auto &a_in = std::get<1>(GetInput());
  const auto &b_in = std::get<2>(GetInput());

  const std::size_t q = ChooseQ(n);
  const std::size_t bs = CeilDiv(n, q);
  const std::size_t np = q * bs;

  std::vector<double> a(np * np, 0.0);
  std::vector<double> b(np * np, 0.0);
  std::vector<double> c(np * np, 0.0);

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      a[Idx(np, i, j)] = a_in[Idx(n, i, j)];
      b[Idx(np, i, j)] = b_in[Idx(n, i, j)];
    }
  }

  for (std::size_t bi = 0; bi < q; ++bi) {
    for (std::size_t bj = 0; bj < q; ++bj) {
      for (std::size_t step = 0; step < q; ++step) {
        const std::size_t bk = (bj + bi + step) % q;
        MulAddBlock(a, b, c, np, bs, bi, bk, bj);
      }
    }
  }

  std::vector<double> out(n * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      out[Idx(n, i, j)] = c[Idx(np, i, j)];
    }
  }

  GetOutput() = std::move(out);
  return true;
}

bool CheremkhinAMatrMultCannonAlgSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace cheremkhin_a_matr_mult_cannon_alg
