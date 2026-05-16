#include "cheremkhin_a_matr_mult_cannon_alg/tbb/include/ops_tbb.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "cheremkhin_a_matr_mult_cannon_alg/common/include/common.hpp"
#include "oneapi/tbb/blocked_range2d.h"
#include "oneapi/tbb/global_control.h"
#include "oneapi/tbb/parallel_for.h"
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

CheremkhinAMatrMultCannonAlgTBB::CheremkhinAMatrMultCannonAlgTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool CheremkhinAMatrMultCannonAlgTBB::ValidationImpl() {
  const std::size_t n = std::get<0>(GetInput());
  const auto &a = std::get<1>(GetInput());
  const auto &b = std::get<2>(GetInput());
  return n > 0 && a.size() == n * n && b.size() == n * n;
}

bool CheremkhinAMatrMultCannonAlgTBB::PreProcessingImpl() {
  GetOutput() = {};
  return true;
}

bool CheremkhinAMatrMultCannonAlgTBB::RunImpl() {
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

  oneapi::tbb::global_control control(oneapi::tbb::global_control::max_allowed_parallelism, requested_threads);

  oneapi::tbb::parallel_for(std::int64_t{0}, n64, [&](std::int64_t i) {
    const auto row = static_cast<std::size_t>(i);
    for (std::size_t j = 0; j < n; ++j) {
      a[Idx(np, row, j)] = a_in[Idx(n, row, j)];
      b[Idx(np, row, j)] = b_in[Idx(n, row, j)];
    }
  });

  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range2d<std::int64_t>(0, q64, 0, q64),
                            [&](const oneapi::tbb::blocked_range2d<std::int64_t> &range) {
    for (std::int64_t bi = range.rows().begin(); bi != range.rows().end(); ++bi) {
      for (std::int64_t bj = range.cols().begin(); bj != range.cols().end(); ++bj) {
        for (std::size_t step = 0; step < q; ++step) {
          const std::size_t bk = (static_cast<std::size_t>(bi) + static_cast<std::size_t>(bj) + step) % q;
          MulAddBlock(a, b, c, np, bs, static_cast<std::size_t>(bi), bk, static_cast<std::size_t>(bj));
        }
      }
    }
  });

  std::vector<double> out(n * n, 0.0);

  oneapi::tbb::parallel_for(std::int64_t{0}, n64, [&](std::int64_t i) {
    const auto row = static_cast<std::size_t>(i);
    for (std::size_t j = 0; j < n; ++j) {
      out[Idx(n, row, j)] = c[Idx(np, row, j)];
    }
  });

  GetOutput() = std::move(out);
  return true;
}

bool CheremkhinAMatrMultCannonAlgTBB::PostProcessingImpl() {
  return true;
}

}  //  namespace cheremkhin_a_matr_mult_cannon_alg
