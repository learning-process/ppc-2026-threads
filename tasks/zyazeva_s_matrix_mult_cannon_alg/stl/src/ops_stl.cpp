#include "zyazeva_s_matrix_mult_cannon_alg/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <thread>
#include <utility>
#include <vector>

#include "zyazeva_s_matrix_mult_cannon_alg/common/include/common.hpp"

namespace {

void ParallelFor(size_t count, const std::function<void(size_t)> &func) {
  size_t threads_count = std::thread::hardware_concurrency();

  if (threads_count == 0) {
    threads_count = 4;
  }

  threads_count = std::min<size_t>(threads_count, count);

  std::vector<std::thread> threads(threads_count);

  size_t block_size = count / threads_count;
  size_t remainder = count % threads_count;

  size_t begin = 0;

  for (size_t th = 0; th < threads_count; ++th) {
    size_t end = begin + block_size + (th < remainder ? 1 : 0);

    threads[th] = std::thread([begin, end, &func]() {
      for (size_t i = begin; i < end; ++i) {
        func(i);
      }
    });

    begin = end;
  }

  for (auto &th : threads) {
    th.join();
  }
}

void MulBlock(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, size_t bs) {
  for (size_t i = 0; i < bs; ++i) {
    for (size_t k = 0; k < bs; ++k) {
      double v = a[(i * bs) + k];

      for (size_t j = 0; j < bs; ++j) {
        c[(i * bs) + j] += v * b[(k * bs) + j];
      }
    }
  }
}

void RegularMultiplication(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c,
                           size_t n) {
  ParallelFor(n, [&](size_t i) {
    for (size_t j = 0; j < n; ++j) {
      double s = 0.0;

      for (size_t k = 0; k < n; ++k) {
        s += a[(i * n) + k] * b[(k * n) + j];
      }

      c[(i * n) + j] = s;
    }
  });
}

void InitializeBlocks(const std::vector<double> &a, const std::vector<double> &b, std::vector<std::vector<double>> &ba,
                      std::vector<std::vector<double>> &bb, size_t g, size_t bs, size_t n) {
  ParallelFor(g * g, [&](size_t id) {
    size_t i = id / g;
    size_t j = id % g;

    ba[id].assign(bs * bs, 0.0);
    bb[id].assign(bs * bs, 0.0);

    for (size_t bi = 0; bi < bs; ++bi) {
      for (size_t bj = 0; bj < bs; ++bj) {
        size_t gi = (i * bs) + bi;
        size_t gj = (j * bs) + bj;
        size_t li = (bi * bs) + bj;

        ba[id][li] = a[(gi * n) + gj];
        bb[id][li] = b[(gi * n) + gj];
      }
    }
  });
}

void AlignBlocks(const std::vector<std::vector<double>> &ba, const std::vector<std::vector<double>> &bb,
                 std::vector<std::vector<double>> &aa, std::vector<std::vector<double>> &ab, size_t g) {
  ParallelFor(g * g, [&](size_t id) {
    size_t i = id / g;
    size_t j = id % g;

    aa[id] = ba[(i * g) + ((j + i) % g)];
    ab[id] = bb[(((i + j) % g) * g) + j];
  });
}

void CannonStep(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &b,
                std::vector<std::vector<double>> &c, size_t g, size_t bs) {
  ParallelFor(g * g, [&](size_t id) { MulBlock(a[id], b[id], c[id], bs); });
}

void Shift(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &b,
           std::vector<std::vector<double>> &na, std::vector<std::vector<double>> &nb, size_t g) {
  ParallelFor(g * g, [&](size_t id) {
    size_t i = id / g;
    size_t j = id % g;

    na[id] = a[(i * g) + ((j + 1) % g)];
    nb[id] = b[(((i + 1) % g) * g) + j];
  });
}

void Assemble(const std::vector<std::vector<double>> &c, std::vector<double> &r, size_t g, size_t bs, size_t n) {
  ParallelFor(g * g, [&](size_t id) {
    size_t i = id / g;
    size_t j = id % g;

    for (size_t bi = 0; bi < bs; ++bi) {
      for (size_t bj = 0; bj < bs; ++bj) {
        size_t gi = (i * bs) + bi;
        size_t gj = (j * bs) + bj;

        r[(gi * n) + gj] = c[id][(bi * bs) + bj];
      }
    }
  });
}

}  // namespace

namespace zyazeva_s_matrix_mult_cannon_alg {

ZyazevaSMatrixMultCannonAlgSTL::ZyazevaSMatrixMultCannonAlgSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ZyazevaSMatrixMultCannonAlgSTL::ValidationImpl() {
  auto sz = std::get<0>(GetInput());

  auto &a = std::get<1>(GetInput());
  auto &b = std::get<2>(GetInput());

  return sz > 0 && a.size() == static_cast<size_t>(sz * sz) && b.size() == static_cast<size_t>(sz * sz);
}

bool ZyazevaSMatrixMultCannonAlgSTL::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool ZyazevaSMatrixMultCannonAlgSTL::RunImpl() {
  size_t n = std::get<0>(GetInput());

  auto &a = std::get<1>(GetInput());
  auto &b = std::get<2>(GetInput());

  std::vector<double> res(n * n, 0.0);

  auto g = static_cast<size_t>(std::sqrt(n));

  if (g <= 1 || g * g != n || n % g != 0) {
    RegularMultiplication(a, b, res, n);

    GetOutput() = res;
    return true;
  }

  size_t bs = n / g;

  std::vector<std::vector<double>> ba(g * g);
  std::vector<std::vector<double>> bb(g * g);

  std::vector<std::vector<double>> aa(g * g);
  std::vector<std::vector<double>> ab(g * g);

  std::vector<std::vector<double>> c(g * g, std::vector<double>(bs * bs, 0.0));

  InitializeBlocks(a, b, ba, bb, g, bs, n);

  AlignBlocks(ba, bb, aa, ab, g);

  for (size_t step = 0; step < g; ++step) {
    CannonStep(aa, ab, c, g, bs);

    if (step < g - 1) {
      std::vector<std::vector<double>> na(g * g);
      std::vector<std::vector<double>> nb(g * g);

      Shift(aa, ab, na, nb, g);

      aa = std::move(na);
      ab = std::move(nb);
    }
  }

  Assemble(c, res, g, bs, n);

  GetOutput() = res;

  return true;
}

bool ZyazevaSMatrixMultCannonAlgSTL::PostProcessingImpl() {
  return true;
}

}  // namespace zyazeva_s_matrix_mult_cannon_alg
