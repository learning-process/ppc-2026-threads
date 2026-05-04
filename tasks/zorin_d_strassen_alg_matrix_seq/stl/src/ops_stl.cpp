#include "zorin_d_strassen_alg_matrix_seq/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <future>
#include <vector>

#include "util/include/util.hpp"
#include "zorin_d_strassen_alg_matrix_seq/common/include/common.hpp"

namespace zorin_d_strassen_alg_matrix_seq {

namespace {

constexpr std::size_t kCutoff = 256;
constexpr std::size_t kBlockSize = 64;

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

void ZeroMatrix(double *dst, std::size_t stride, std::size_t n) {
  for (std::size_t i = 0; i < n; ++i) {
    std::fill_n(dst + (i * stride), n, 0.0);
  }
}

void AddToBuffer(const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *dst,
                 std::size_t n, double b_coeff) {
  for (std::size_t i = 0; i < n; ++i) {
    const double *a_row = a + (i * a_stride);
    const double *b_row = b + (i * b_stride);
    double *dst_row = dst + (i * n);
    for (std::size_t j = 0; j < n; ++j) {
      dst_row[j] = a_row[j] + (b_coeff * b_row[j]);
    }
  }
}

void MulMicroBlock(const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *c,
                   std::size_t c_stride, std::size_t i_begin, std::size_t i_end, std::size_t k_begin, std::size_t k_end,
                   std::size_t j_begin, std::size_t j_end) {
  for (std::size_t i = i_begin; i < i_end; ++i) {
    double *c_row = c + (i * c_stride);
    const double *a_row = a + (i * a_stride);
    for (std::size_t k = k_begin; k < k_end; ++k) {
      const double aik = a_row[k];
      const double *b_row = b + (k * b_stride);
      for (std::size_t j = j_begin; j < j_end; ++j) {
        c_row[j] += aik * b_row[j];
      }
    }
  }
}

void MulJBlocks(const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *c,
                std::size_t c_stride, std::size_t n, std::size_t ii, std::size_t i_end, std::size_t kk,
                std::size_t k_end) {
  for (std::size_t jj = 0; jj < n; jj += kBlockSize) {
    const std::size_t j_end = std::min(jj + kBlockSize, n);
    MulMicroBlock(a, a_stride, b, b_stride, c, c_stride, ii, i_end, kk, k_end, jj, j_end);
  }
}

void MulKBlocks(const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *c,
                std::size_t c_stride, std::size_t n, std::size_t ii, std::size_t i_end) {
  for (std::size_t kk = 0; kk < n; kk += kBlockSize) {
    const std::size_t k_end = std::min(kk + kBlockSize, n);
    MulJBlocks(a, a_stride, b, b_stride, c, c_stride, n, ii, i_end, kk, k_end);
  }
}

void NaiveMulBlocked(const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *c,
                     std::size_t c_stride, std::size_t n) {
  ZeroMatrix(c, c_stride, n);

  for (std::size_t ii = 0; ii < n; ii += kBlockSize) {
    const std::size_t i_end = std::min(ii + kBlockSize, n);
    MulKBlocks(a, a_stride, b, b_stride, c, c_stride, n, ii, i_end);
  }
}

void CombineQuadrants(const std::vector<double> &m1, const std::vector<double> &m2, const std::vector<double> &m3,
                      const std::vector<double> &m4, const std::vector<double> &m5, const std::vector<double> &m6,
                      const std::vector<double> &m7, double *c, std::size_t c_stride, std::size_t half) {
  for (std::size_t i = 0; i < half; ++i) {
    double *c11 = c + (i * c_stride);
    double *c12 = c11 + half;
    double *c21 = c + ((i + half) * c_stride);
    double *c22 = c21 + half;

    const double *m1_row = m1.data() + (i * half);
    const double *m2_row = m2.data() + (i * half);
    const double *m3_row = m3.data() + (i * half);
    const double *m4_row = m4.data() + (i * half);
    const double *m5_row = m5.data() + (i * half);
    const double *m6_row = m6.data() + (i * half);
    const double *m7_row = m7.data() + (i * half);

    for (std::size_t j = 0; j < half; ++j) {
      c11[j] = m1_row[j] + m4_row[j] - m5_row[j] + m7_row[j];
      c12[j] = m3_row[j] + m5_row[j];
      c21[j] = m2_row[j] + m4_row[j];
      c22[j] = m1_row[j] - m2_row[j] + m3_row[j] + m6_row[j];
    }
  }
}

using StrassenSeqFn = void (*)(const double *, std::size_t, const double *, std::size_t, double *, std::size_t,
                               std::size_t);

void StrassenSeqImpl(const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *c,
                     std::size_t c_stride, std::size_t n);

constexpr StrassenSeqFn kStrassenSeqFn = &StrassenSeqImpl;

void ComputeProduct(const double *a1, std::size_t a1_stride, const double *a2, std::size_t a2_stride, double a2_coeff,
                    const double *b1, std::size_t b1_stride, const double *b2, std::size_t b2_stride, double b2_coeff,
                    std::vector<double> &out, std::size_t n) {
  std::vector<double> lhs(n * n);
  std::vector<double> rhs(n * n);
  AddToBuffer(a1, a1_stride, a2, a2_stride, lhs.data(), n, a2_coeff);
  AddToBuffer(b1, b1_stride, b2, b2_stride, rhs.data(), n, b2_coeff);
  out.assign(n * n, 0.0);
  kStrassenSeqFn(lhs.data(), n, rhs.data(), n, out.data(), n, n);
}

void ComputeProductSingle(const double *a, std::size_t a_stride, const double *b1, std::size_t b1_stride,
                          const double *b2, std::size_t b2_stride, double b2_coeff, std::vector<double> &out,
                          std::size_t n) {
  std::vector<double> rhs(n * n);
  AddToBuffer(b1, b1_stride, b2, b2_stride, rhs.data(), n, b2_coeff);
  out.assign(n * n, 0.0);
  kStrassenSeqFn(a, a_stride, rhs.data(), n, out.data(), n, n);
}

void ComputeProductSingleLeft(const double *a1, std::size_t a1_stride, const double *a2, std::size_t a2_stride,
                              double a2_coeff, const double *b, std::size_t b_stride, std::vector<double> &out,
                              std::size_t n) {
  std::vector<double> lhs(n * n);
  AddToBuffer(a1, a1_stride, a2, a2_stride, lhs.data(), n, a2_coeff);
  out.assign(n * n, 0.0);
  kStrassenSeqFn(lhs.data(), n, b, b_stride, out.data(), n, n);
}

void StrassenTopStl(const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *c,
                    std::size_t c_stride, std::size_t n) {
  if (n <= kCutoff || ppc::util::GetNumThreads() <= 1) {
    kStrassenSeqFn(a, a_stride, b, b_stride, c, c_stride, n);
    return;
  }

  const std::size_t half = n / 2;

  const double *a11 = a;
  const double *a12 = a + half;
  const double *a21 = a + (half * a_stride);
  const double *a22 = a21 + half;

  const double *b11 = b;
  const double *b12 = b + half;
  const double *b21 = b + (half * b_stride);
  const double *b22 = b21 + half;

  std::vector<double> m1;
  std::vector<double> m2;
  std::vector<double> m3;
  std::vector<double> m4;
  std::vector<double> m5;
  std::vector<double> m6;
  std::vector<double> m7;

  auto f1 = std::async(std::launch::async, [&] {
    ComputeProduct(a11, a_stride, a22, a_stride, 1.0, b11, b_stride, b22, b_stride, 1.0, m1, half);
  });
  auto f2 = std::async(std::launch::async,
                       [&] { ComputeProductSingleLeft(a21, a_stride, a22, a_stride, 1.0, b11, b_stride, m2, half); });
  auto f3 = std::async(std::launch::async,
                       [&] { ComputeProductSingle(a11, a_stride, b12, b_stride, b22, b_stride, -1.0, m3, half); });
  auto f4 = std::async(std::launch::async,
                       [&] { ComputeProductSingle(a22, a_stride, b21, b_stride, b11, b_stride, -1.0, m4, half); });
  auto f5 = std::async(std::launch::async,
                       [&] { ComputeProductSingleLeft(a11, a_stride, a12, a_stride, 1.0, b22, b_stride, m5, half); });
  auto f6 = std::async(std::launch::async, [&] {
    ComputeProduct(a21, a_stride, a11, a_stride, -1.0, b11, b_stride, b12, b_stride, 1.0, m6, half);
  });
  auto f7 = std::async(std::launch::async, [&] {
    ComputeProduct(a12, a_stride, a22, a_stride, -1.0, b21, b_stride, b22, b_stride, 1.0, m7, half);
  });

  f1.get();
  f2.get();
  f3.get();
  f4.get();
  f5.get();
  f6.get();
  f7.get();

  CombineQuadrants(m1, m2, m3, m4, m5, m6, m7, c, c_stride, half);
}

void StrassenSeqImpl(const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *c,
                     std::size_t c_stride, std::size_t n) {
  if (n <= kCutoff) {
    NaiveMulBlocked(a, a_stride, b, b_stride, c, c_stride, n);
    return;
  }

  const std::size_t half = n / 2;

  const double *a11 = a;
  const double *a12 = a + half;
  const double *a21 = a + (half * a_stride);
  const double *a22 = a21 + half;

  const double *b11 = b;
  const double *b12 = b + half;
  const double *b21 = b + (half * b_stride);
  const double *b22 = b21 + half;

  std::vector<double> m1(half * half);
  std::vector<double> m2(half * half);
  std::vector<double> m3(half * half);
  std::vector<double> m4(half * half);
  std::vector<double> m5(half * half);
  std::vector<double> m6(half * half);
  std::vector<double> m7(half * half);
  std::vector<double> lhs(half * half);
  std::vector<double> rhs(half * half);

  AddToBuffer(a11, a_stride, a22, a_stride, lhs.data(), half, 1.0);
  AddToBuffer(b11, b_stride, b22, b_stride, rhs.data(), half, 1.0);
  kStrassenSeqFn(lhs.data(), half, rhs.data(), half, m1.data(), half, half);

  AddToBuffer(a21, a_stride, a22, a_stride, lhs.data(), half, 1.0);
  kStrassenSeqFn(lhs.data(), half, b11, b_stride, m2.data(), half, half);

  AddToBuffer(b12, b_stride, b22, b_stride, rhs.data(), half, -1.0);
  kStrassenSeqFn(a11, a_stride, rhs.data(), half, m3.data(), half, half);

  AddToBuffer(b21, b_stride, b11, b_stride, rhs.data(), half, -1.0);
  kStrassenSeqFn(a22, a_stride, rhs.data(), half, m4.data(), half, half);

  AddToBuffer(a11, a_stride, a12, a_stride, lhs.data(), half, 1.0);
  kStrassenSeqFn(lhs.data(), half, b22, b_stride, m5.data(), half, half);

  AddToBuffer(a21, a_stride, a11, a_stride, lhs.data(), half, -1.0);
  AddToBuffer(b11, b_stride, b12, b_stride, rhs.data(), half, 1.0);
  kStrassenSeqFn(lhs.data(), half, rhs.data(), half, m6.data(), half, half);

  AddToBuffer(a12, a_stride, a22, a_stride, lhs.data(), half, -1.0);
  AddToBuffer(b21, b_stride, b22, b_stride, rhs.data(), half, 1.0);
  kStrassenSeqFn(lhs.data(), half, rhs.data(), half, m7.data(), half, half);

  CombineQuadrants(m1, m2, m3, m4, m5, m6, m7, c, c_stride, half);
}

}  // namespace

ZorinDStrassenAlgMatrixSTL::ZorinDStrassenAlgMatrixSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ZorinDStrassenAlgMatrixSTL::ValidationImpl() {
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

bool ZorinDStrassenAlgMatrixSTL::PreProcessingImpl() {
  const auto n = GetInput().n;
  GetOutput().assign(n * n, 0.0);
  return true;
}

bool ZorinDStrassenAlgMatrixSTL::RunImpl() {
  const auto &in = GetInput();
  const std::size_t n = in.n;
  const std::size_t padded = NextPow2(n);

  std::vector<double> a_pad(padded * padded, 0.0);
  std::vector<double> b_pad(padded * padded, 0.0);
  std::vector<double> c_pad(padded * padded, 0.0);

  for (std::size_t i = 0; i < n; ++i) {
    std::copy_n(in.a.data() + (i * n), n, a_pad.data() + (i * padded));
    std::copy_n(in.b.data() + (i * n), n, b_pad.data() + (i * padded));
  }

  StrassenTopStl(a_pad.data(), padded, b_pad.data(), padded, c_pad.data(), padded, padded);

  auto &out = GetOutput();
  for (std::size_t i = 0; i < n; ++i) {
    std::copy_n(c_pad.data() + (i * padded), n, out.data() + (i * n));
  }

  return true;
}

bool ZorinDStrassenAlgMatrixSTL::PostProcessingImpl() {
  return true;
}

}  // namespace zorin_d_strassen_alg_matrix_seq
