#include "zyazeva_s_matrix_mult_cannon_alg/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <execution>
#include <thread>
#include <vector>

namespace {

std::vector<int> MakeIdx(int n) {
  std::vector<int> v(n);
  for (int i = 0; i < n; ++i) {
    v[i] = i;
  }
  return v;
}

void MulBlock(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, int bs) {
  for (int i = 0; i < bs; ++i) {
    for (int k = 0; k < bs; ++k) {
      double v = a[i * bs + k];
      for (int j = 0; j < bs; ++j) {
        c[i * bs + j] += v * b[k * bs + j];
      }
    }
  }
}

void RegularMultiplication(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, int n) {
  auto idx = MakeIdx(n);

  std::for_each(std::execution::par, idx.begin(), idx.end(), [&](int i) {
    for (int j = 0; j < n; ++j) {
      double s = 0;
      for (int k = 0; k < n; ++k) {
        s += a[i * n + k] * b[k * n + j];
      }
      c[i * n + j] = s;
    }
  });
}

void InitializeBlocks(const std::vector<double> &a, const std::vector<double> &b, std::vector<std::vector<double>> &ba,
                      std::vector<std::vector<double>> &bb, int g, int bs, int n) {
  for (int i = 0; i < g; ++i) {
    for (int j = 0; j < g; ++j) {
      int id = i * g + j;
      ba[id].assign(bs * bs, 0.0);
      bb[id].assign(bs * bs, 0.0);

      for (int bi = 0; bi < bs; ++bi) {
        for (int bj = 0; bj < bs; ++bj) {
          int gi = i * bs + bi;
          int gj = j * bs + bj;
          int li = bi * bs + bj;

          ba[id][li] = a[gi * n + gj];
          bb[id][li] = b[gi * n + gj];
        }
      }
    }
  }
}

void AlignBlocks(const std::vector<std::vector<double>> &ba, const std::vector<std::vector<double>> &bb,
                 std::vector<std::vector<double>> &aa, std::vector<std::vector<double>> &ab, int g) {
  auto idx = MakeIdx(g * g);

  std::for_each(std::execution::par, idx.begin(), idx.end(), [&](int id) {
    int i = id / g, j = id % g;
    aa[id] = ba[i * g + (j + i) % g];
    ab[id] = bb[((i + j) % g) * g + j];
  });
}

void CannonStep(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &b,
                std::vector<std::vector<double>> &c, int g, int bs) {
  auto idx = MakeIdx(g * g);

  std::for_each(std::execution::par, idx.begin(), idx.end(), [&](int id) { MulBlock(a[id], b[id], c[id], bs); });
}

void Shift(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &b,
           std::vector<std::vector<double>> &na, std::vector<std::vector<double>> &nb, int g) {
  auto idx = MakeIdx(g * g);

  std::for_each(std::execution::par, idx.begin(), idx.end(), [&](int id) {
    int i = id / g, j = id % g;
    na[id] = a[i * g + (j + 1) % g];
    nb[id] = b[((i + 1) % g) * g + j];
  });
}

void Assemble(const std::vector<std::vector<double>> &c, std::vector<double> &r, int g, int bs, int n) {
  auto idx = MakeIdx(g * g);

  std::for_each(std::execution::par, idx.begin(), idx.end(), [&](int id) {
    int i = id / g, j = id % g;

    for (int bi = 0; bi < bs; ++bi) {
      for (int bj = 0; bj < bs; ++bj) {
        int gi = i * bs + bi;
        int gj = j * bs + bj;
        r[gi * n + gj] = c[id][bi * bs + bj];
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
  return sz > 0 && a.size() == sz * sz && b.size() == sz * sz;
}

bool ZyazevaSMatrixMultCannonAlgSTL::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool ZyazevaSMatrixMultCannonAlgSTL::RunImpl() {
  int n = std::get<0>(GetInput());
  auto &a = std::get<1>(GetInput());
  auto &b = std::get<2>(GetInput());

  std::vector<double> res(n * n, 0.0);

  int g = static_cast<int>(std::sqrt(n));

  if (g <= 1 || g * g != n || n % g != 0) {
    RegularMultiplication(a, b, res, n);
    GetOutput() = res;
    return true;
  }

  int bs = n / g;

  std::vector<std::vector<double>> ba(g * g), bb(g * g);
  std::vector<std::vector<double>> aa(g * g), ab(g * g);
  std::vector<std::vector<double>> c(g * g, std::vector<double>(bs * bs, 0.0));

  InitializeBlocks(a, b, ba, bb, g, bs, n);
  AlignBlocks(ba, bb, aa, ab, g);

  for (int i = 0; i < g; ++i) {
    CannonStep(aa, ab, c, g, bs);

    if (i < g - 1) {
      std::vector<std::vector<double>> na(g * g), nb(g * g);
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
