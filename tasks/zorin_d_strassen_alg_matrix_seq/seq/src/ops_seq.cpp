#include "zorin_d_strassen_alg_matrix_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

#include "zorin_d_strassen_alg_matrix_seq/common/include/common.hpp"

namespace zorin_d_strassen_alg_matrix_seq {

namespace {

constexpr std::size_t kCutoff = 64;

std::size_t NextPow2(std::size_t x) {
  if (x <= 1) {
    return 1;
  }
  std::size_t p = 1;
  while (p < x) {
    p <<= 1;
  }
  return p;
}

void NaiveMul(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, std::size_t n) {
  std::fill(c.begin(), c.end(), 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t i_row = i * n;
    for (std::size_t k = 0; k < n; ++k) {
      const double aik = a[i_row + k];
      const std::size_t k_row = k * n;
      for (std::size_t j = 0; j < n; ++j) {
        c[i_row + j] += aik * b[k_row + j];
      }
    }
  }
}

std::vector<double> AddVec(const std::vector<double> &x, const std::vector<double> &y) {
  assert(x.size() == y.size());
  std::vector<double> r(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    r[i] = x[i] + y[i];
  }
  return r;
}

std::vector<double> SubVec(const std::vector<double> &x, const std::vector<double> &y) {
  assert(x.size() == y.size());
  std::vector<double> r(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    r[i] = x[i] - y[i];
  }
  return r;
}

void SplitMat(const std::vector<double> &m, std::size_t n, std::vector<double> &m11, std::vector<double> &m12,
              std::vector<double> &m21, std::vector<double> &m22) {
  const std::size_t k = n / 2;
  m11.assign(k * k, 0.0);
  m12.assign(k * k, 0.0);
  m21.assign(k * k, 0.0);
  m22.assign(k * k, 0.0);

  for (std::size_t i = 0; i < k; ++i) {
    for (std::size_t j = 0; j < k; ++j) {
      m11[(i * k) + j] = m[(i * n) + j];
      m12[(i * k) + j] = m[(i * n) + (j + k)];
      m21[(i * k) + j] = m[((i + k) * n) + j];
      m22[(i * k) + j] = m[((i + k) * n) + (j + k)];
    }
  }
}

std::vector<double> JoinMat(const std::vector<double> &c11, const std::vector<double> &c12,
                            const std::vector<double> &c21, const std::vector<double> &c22, std::size_t n) {
  const std::size_t k = n / 2;
  std::vector<double> c(n * n, 0.0);
  for (std::size_t i = 0; i < k; ++i) {
    for (std::size_t j = 0; j < k; ++j) {
      c[(i * n) + j] = c11[(i * k) + j];
      c[(i * n) + (j + k)] = c12[(i * k) + j];
      c[((i + k) * n) + j] = c21[(i * k) + j];
      c[((i + k) * n) + (j + k)] = c22[(i * k) + j];
    }
  }
  return c;
}

struct Frame {
  std::size_t n{};
  std::size_t half{};
  int stage{};
  int parent_slot{-1};

  std::vector<double> a;
  std::vector<double> b;
  std::vector<double> c;

  std::vector<double> a11;
  std::vector<double> a12;
  std::vector<double> a21;
  std::vector<double> a22;

  std::vector<double> b11;
  std::vector<double> b12;
  std::vector<double> b21;
  std::vector<double> b22;

  std::vector<double> m1;
  std::vector<double> m2;
  std::vector<double> m3;
  std::vector<double> m4;
  std::vector<double> m5;
  std::vector<double> m6;
  std::vector<double> m7;

  Frame(std::vector<double> a_in, std::vector<double> b_in, std::size_t n_in, int parent_slot_in)
      : n(n_in),
        half(n_in / 2),
        stage(0),
        parent_slot(parent_slot_in),
        a(std::move(a_in)),
        b(std::move(b_in)),
        c(n_in * n_in, 0.0) {}
};

std::vector<double> &SelectM(Frame &f, int slot) {
  if (slot == 1) {
    return f.m1;
  }
  if (slot == 2) {
    return f.m2;
  }
  if (slot == 3) {
    return f.m3;
  }
  if (slot == 4) {
    return f.m4;
  }
  if (slot == 5) {
    return f.m5;
  }
  if (slot == 6) {
    return f.m6;
  }
  return f.m7;
}

Frame MakeChildForSlot(const Frame &f, int slot) {
  const std::size_t k = f.half;
  if (slot == 1) {
    return {AddVec(f.a11, f.a22), AddVec(f.b11, f.b22), k, 1};
  }
  if (slot == 2) {
    return {AddVec(f.a21, f.a22), f.b11, k, 2};
  }
  if (slot == 3) {
    return {f.a11, SubVec(f.b12, f.b22), k, 3};
  }
  if (slot == 4) {
    return {f.a22, SubVec(f.b21, f.b11), k, 4};
  }
  if (slot == 5) {
    return {AddVec(f.a11, f.a12), f.b22, k, 5};
  }
  if (slot == 6) {
    return {SubVec(f.a21, f.a11), AddVec(f.b11, f.b12), k, 6};
  }
  return {SubVec(f.a12, f.a22), AddVec(f.b21, f.b22), k, 7};
}

std::vector<double> StrassenIter(std::vector<double> a, std::vector<double> b, std::size_t n) {
  std::vector<Frame> st;
  st.emplace_back(std::move(a), std::move(b), n, -1);

  std::vector<double> root_out;

  while (!st.empty()) {
    Frame &f = st.back();

    if (f.n <= kCutoff) {
      NaiveMul(f.a, f.b, f.c, f.n);

      Frame finished = std::move(f);
      st.pop_back();

      if (st.empty()) {
        root_out = std::move(finished.c);
        continue;
      }

      Frame &parent = st.back();
      SelectM(parent, finished.parent_slot) = std::move(finished.c);
      parent.stage += 1;
      continue;
    }

    if (f.stage == 0) {
      SplitMat(f.a, f.n, f.a11, f.a12, f.a21, f.a22);
      SplitMat(f.b, f.n, f.b11, f.b12, f.b21, f.b22);
      st.push_back(MakeChildForSlot(f, 1));
      continue;
    }

    if (f.stage >= 1 && f.stage <= 6) {
      const int next_slot = f.stage + 1;
      st.push_back(MakeChildForSlot(f, next_slot));
      continue;
    }

    const auto c11 = AddVec(SubVec(AddVec(f.m1, f.m4), f.m5), f.m7);
    const auto c12 = AddVec(f.m3, f.m5);
    const auto c21 = AddVec(f.m2, f.m4);
    const auto c22 = AddVec(AddVec(SubVec(f.m1, f.m2), f.m3), f.m6);

    f.c = JoinMat(c11, c12, c21, c22, f.n);

    Frame finished = std::move(f);
    st.pop_back();

    if (st.empty()) {
      root_out = std::move(finished.c);
      continue;
    }

    Frame &parent = st.back();
    SelectM(parent, finished.parent_slot) = std::move(finished.c);
    parent.stage += 1;
  }

  return root_out;
}

}  // namespace

ZorinDStrassenAlgMatrixSEQ::ZorinDStrassenAlgMatrixSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ZorinDStrassenAlgMatrixSEQ::ValidationImpl() {
  const auto &in = GetInput();
  if (in.n == 0) {
    return false;
  }
  if (in.a.size() != in.n * in.n) {
    return false;
  }
  if (in.b.size() != in.n * in.n) {
    return false;
  }
  if (!GetOutput().empty()) {
    return false;
  }
  return true;
}

bool ZorinDStrassenAlgMatrixSEQ::PreProcessingImpl() {
  const auto n = GetInput().n;
  GetOutput().assign(n * n, 0.0);
  return true;
}

bool ZorinDStrassenAlgMatrixSEQ::RunImpl() {
  const auto &in = GetInput();
  const std::size_t n = in.n;
  const std::size_t m = NextPow2(n);

  std::vector<double> a_pad(m * m, 0.0);
  std::vector<double> b_pad(m * m, 0.0);

  for (std::size_t i = 0; i < n; ++i) {
    std::copy_n(&in.a[i * n], n, &a_pad[i * m]);
    std::copy_n(&in.b[i * n], n, &b_pad[i * m]);
  }

  const auto c_pad = StrassenIter(std::move(a_pad), std::move(b_pad), m);

  auto &out = GetOutput();
  for (std::size_t i = 0; i < n; ++i) {
    std::copy_n(&c_pad[i * m], n, &out[i * n]);
  }

  return true;
}

bool ZorinDStrassenAlgMatrixSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace zorin_d_strassen_alg_matrix_seq
