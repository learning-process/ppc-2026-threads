#include "zorin_d_strassen_alg_matrix_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>

namespace zorin_d_strassen_alg_matrix_seq {

namespace {

constexpr std::size_t k_cutoff = 64;

inline std::size_t next_pow2(std::size_t x) {
  if (x <= 1) {
    return 1;
  }
  std::size_t p = 1;
  while (p < x) {
    p <<= 1;
  }
  return p;
}

inline void naive_mul(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c,
                      std::size_t n) {
  std::fill(c.begin(), c.end(), 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t k = 0; k < n; ++k) {
      const double aik = a[(i * n) + k];
      for (std::size_t j = 0; j < n; ++j) {
        c[(i * n) + j] += aik * b[(k * n) + j];
      }
    }
  }
}

inline std::vector<double> add_vec(const std::vector<double> &x, const std::vector<double> &y) {
  assert(x.size() == y.size());
  std::vector<double> r(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    r[i] = x[i] + y[i];
  }
  return r;
}

inline std::vector<double> sub_vec(const std::vector<double> &x, const std::vector<double> &y) {
  assert(x.size() == y.size());
  std::vector<double> r(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    r[i] = x[i] - y[i];
  }
  return r;
}

inline void split_mat(const std::vector<double> &m, std::size_t n, std::vector<double> &m11, std::vector<double> &m12,
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

inline std::vector<double> join_mat(const std::vector<double> &c11, const std::vector<double> &c12,
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

std::vector<double> strassen_rec(const std::vector<double> &a, const std::vector<double> &b, std::size_t n) {
  if (n <= k_cutoff) {
    std::vector<double> c(n * n);
    naive_mul(a, b, c, n);
    return c;
  }

  const std::size_t k = n / 2;

  std::vector<double> a11;
  std::vector<double> a12;
  std::vector<double> a21;
  std::vector<double> a22;

  std::vector<double> b11;
  std::vector<double> b12;
  std::vector<double> b21;
  std::vector<double> b22;

  split_mat(a, n, a11, a12, a21, a22);
  split_mat(b, n, b11, b12, b21, b22);

  const auto m1 = strassen_rec(add_vec(a11, a22), add_vec(b11, b22), k);
  const auto m2 = strassen_rec(add_vec(a21, a22), b11, k);
  const auto m3 = strassen_rec(a11, sub_vec(b12, b22), k);
  const auto m4 = strassen_rec(a22, sub_vec(b21, b11), k);
  const auto m5 = strassen_rec(add_vec(a11, a12), b22, k);
  const auto m6 = strassen_rec(sub_vec(a21, a11), add_vec(b11, b12), k);
  const auto m7 = strassen_rec(sub_vec(a12, a22), add_vec(b21, b22), k);

  const auto c11 = add_vec(sub_vec(add_vec(m1, m4), m5), m7);
  const auto c12 = add_vec(m3, m5);
  const auto c21 = add_vec(m2, m4);
  const auto c22 = add_vec(add_vec(sub_vec(m1, m2), m3), m6);

  return join_mat(c11, c12, c21, c22, n);
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
  if (in.A.size() != in.n * in.n) {
    return false;
  }
  if (in.B.size() != in.n * in.n) {
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
  const std::size_t m = next_pow2(n);

  std::vector<double> a_pad(m * m, 0.0);
  std::vector<double> b_pad(m * m, 0.0);

  for (std::size_t i = 0; i < n; ++i) {
    std::copy_n(&in.A[i * n], n, &a_pad[i * m]);
    std::copy_n(&in.B[i * n], n, &b_pad[i * m]);
  }

  const auto c_pad = strassen_rec(a_pad, b_pad, m);

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
