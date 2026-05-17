#include "cheremkhin_a_matr_mult_cannon_alg/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <future>
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

  for (std::size_t ii = 0; ii < bs; ++ii) {
    const std::size_t i = i0 + ii;
    const std::size_t a_row = i * n;
    const std::size_t c_row = i * n;
    double *c_block = c.data() + c_row + j0;

    for (std::size_t kk = 0; kk < bs; ++kk) {
      const std::size_t k = k0 + kk;
      const double aik = a[a_row + k];
      const double *b_block = b.data() + (k * n) + j0;
      for (std::size_t jj = 0; jj < bs; ++jj) {
        c_block[jj] += aik * b_block[jj];
      }
    }
  }
}

template <class Func>
void ParallelFor(std::size_t count, std::size_t requested_threads, Func fn) {
  if (count == 0) {
    return;
  }

  const std::size_t workers = std::max<std::size_t>(1, std::min(count, requested_threads));
  if (workers == 1) {
    for (std::size_t idx = 0; idx < count; ++idx) {
      fn(idx);
    }
    return;
  }

  const std::size_t chunk = CeilDiv(count, workers);
  std::vector<std::future<void>> tasks;
  tasks.reserve(workers);

  for (std::size_t worker = 0; worker < workers; ++worker) {
    const std::size_t begin = worker * chunk;
    const std::size_t end = std::min(begin + chunk, count);
    if (begin >= end) {
      break;
    }

    tasks.emplace_back(std::async(std::launch::async, [begin, end, &fn] {
      for (std::size_t idx = begin; idx < end; ++idx) {
        fn(idx);
      }
    }));
  }

  for (auto &task : tasks) {
    task.get();
  }
}

}  // namespace

CheremkhinAMatrMultCannonAlgSTL::CheremkhinAMatrMultCannonAlgSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool CheremkhinAMatrMultCannonAlgSTL::ValidationImpl() {
  const std::size_t n = std::get<0>(GetInput());
  const auto &a = std::get<1>(GetInput());
  const auto &b = std::get<2>(GetInput());
  return n > 0 && a.size() == n * n && b.size() == n * n;
}

bool CheremkhinAMatrMultCannonAlgSTL::PreProcessingImpl() {
  GetOutput() = {};
  return true;
}

bool CheremkhinAMatrMultCannonAlgSTL::RunImpl() {
  const std::size_t n = std::get<0>(GetInput());
  const auto &a_in = std::get<1>(GetInput());
  const auto &b_in = std::get<2>(GetInput());
  const std::size_t threads = static_cast<std::size_t>(std::max(1, ppc::util::GetNumThreads()));

  const std::size_t q = ChooseQ(n);
  const std::size_t bs = CeilDiv(n, q);
  const std::size_t np = q * bs;

  std::vector<double> a(np * np, 0.0);
  std::vector<double> b(np * np, 0.0);
  std::vector<double> c(np * np, 0.0);

  ParallelFor(n, threads, [&](std::size_t i) {
    for (std::size_t j = 0; j < n; ++j) {
      a[Idx(np, i, j)] = a_in[Idx(n, i, j)];
      b[Idx(np, i, j)] = b_in[Idx(n, i, j)];
    }
  });

  ParallelFor(q * q, threads, [&](std::size_t block_idx) {
    const std::size_t bi = block_idx / q;
    const std::size_t bj = block_idx % q;

    for (std::size_t step = 0; step < q; ++step) {
      const std::size_t bk = (bi + bj + step) % q;
      MulAddBlock(a, b, c, np, bs, bi, bk, bj);
    }
  });

  std::vector<double> out(n * n, 0.0);
  ParallelFor(n, threads, [&](std::size_t i) {
    for (std::size_t j = 0; j < n; ++j) {
      out[Idx(n, i, j)] = c[Idx(np, i, j)];
    }
  });

  GetOutput() = std::move(out);
  return true;
}

bool CheremkhinAMatrMultCannonAlgSTL::PostProcessingImpl() {
  return true;
}

}  // namespace  cheremkhin_a_matr_mult_cannon_alg
