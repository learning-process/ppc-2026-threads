#include "cheremkhin_a_matr_mult_cannon_alg/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "cheremkhin_a_matr_mult_cannon_alg/common/include/common.hpp"
#include "util/include/util.hpp"

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
  const auto bs64 = static_cast<std::int64_t>(bs);

  for (std::size_t ii = 0; ii < bs; ++ii) {
    const std::size_t i = i0 + ii;
    const std::size_t a_row = i * n;
    const std::size_t c_row = i * n;
    double *c_block = c.data() + c_row + j0;
    for (std::size_t kk = 0; kk < bs; ++kk) {
      const std::size_t k = k0 + kk;
      const double aik = a[a_row + k];
      const double *b_block = b.data() + (k * n) + j0;
      for (std::int64_t jj = 0; jj < bs64; ++jj) {
        c_block[jj] += aik * b_block[jj];
      }
    }
  }
}

}  // namespace

CheremkhinAMatrMultCannonAlgOMP::CheremkhinAMatrMultCannonAlgOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool CheremkhinAMatrMultCannonAlgOMP::ValidationImpl() {
  const std::size_t n = std::get<0>(GetInput());
  const auto &a = std::get<1>(GetInput());
  const auto &b = std::get<2>(GetInput());
  return n > 0 && a.size() == n * n && b.size() == n * n;
}

bool CheremkhinAMatrMultCannonAlgOMP::PreProcessingImpl() {
  GetOutput() = {};
  return true;
}

bool CheremkhinAMatrMultCannonAlgOMP::RunImpl() {
  const std::size_t n = std::get<0>(GetInput());
  const auto &a_in = std::get<1>(GetInput());
  const auto &b_in = std::get<2>(GetInput());
  const int requested_threads = ppc::util::GetNumThreads();

  const std::size_t q = ChooseQ(n);
  const std::size_t bs = CeilDiv(n, q);
  const std::size_t np = q * bs;
  const auto n64 = static_cast<std::int64_t>(n);
  const auto q64 = static_cast<std::int64_t>(q);

  std::vector<double> a(np * np, 0.0);
  std::vector<double> b(np * np, 0.0);
  std::vector<double> c(np * np, 0.0);

  omp_set_num_threads(requested_threads);

#pragma omp parallel for default(none) schedule(static) shared(a, b, a_in, b_in, n, np, n64)
  for (std::int64_t i = 0; i < n64; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      a[Idx(np, static_cast<std::size_t>(i), j)] = a_in[Idx(n, static_cast<std::size_t>(i), j)];
      b[Idx(np, static_cast<std::size_t>(i), j)] = b_in[Idx(n, static_cast<std::size_t>(i), j)];
    }
  }

#pragma omp parallel for default(none) collapse(2) schedule(static) shared(a, b, c, np, bs, q, q64)
  for (std::int64_t bi = 0; bi < q64; ++bi) {
    for (std::int64_t bj = 0; bj < q64; ++bj) {
      for (std::size_t step = 0; step < q; ++step) {
        const std::size_t bk = (static_cast<std::size_t>(bi) + static_cast<std::size_t>(bj) + step) % q;
        MulAddBlock(a, b, c, np, bs, static_cast<std::size_t>(bi), bk, static_cast<std::size_t>(bj));
      }
    }
  }

  std::vector<double> out(n * n, 0.0);

#pragma omp parallel for default(none) schedule(static) shared(out, c, n, np, n64)
  for (std::int64_t i = 0; i < n64; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      out[Idx(n, static_cast<std::size_t>(i), j)] = c[Idx(np, static_cast<std::size_t>(i), j)];
    }
  }

  GetOutput() = std::move(out);
  return true;
}

bool CheremkhinAMatrMultCannonAlgOMP::PostProcessingImpl() {
  return true;
}

}  // namespace  cheremkhin_a_matr_mult_cannon_alg
