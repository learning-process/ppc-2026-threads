#include "muhammadkhon_i_stressen_alg/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <future>
#include <vector>

#include "muhammadkhon_i_stressen_alg/common/include/common.hpp"
#include "util/include/util.hpp"

namespace muhammadkhon_i_stressen_alg {

namespace {

constexpr std::size_t kCutoff = 64;
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

void MulKBlocks(const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *c,
                std::size_t c_stride, std::size_t n, std::size_t ii, std::size_t i_end) {
  for (std::size_t kk = 0; kk < n; kk += kBlockSize) {
    const std::size_t k_end = std::min(kk + kBlockSize, n);
    for (std::size_t jj = 0; jj < n; jj += kBlockSize) {
      const std::size_t j_end = std::min(jj + kBlockSize, n);
      MulMicroBlock(a, a_stride, b, b_stride, c, c_stride, ii, i_end, kk, k_end, jj, j_end);
    }
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
    const double *m1r = m1.data() + (i * half);
    const double *m2r = m2.data() + (i * half);
    const double *m3r = m3.data() + (i * half);
    const double *m4r = m4.data() + (i * half);
    const double *m5r = m5.data() + (i * half);
    const double *m6r = m6.data() + (i * half);
    const double *m7r = m7.data() + (i * half);
    for (std::size_t j = 0; j < half; ++j) {
      c11[j] = m1r[j] + m4r[j] - m5r[j] + m7r[j];
      c12[j] = m3r[j] + m5r[j];
      c21[j] = m2r[j] + m4r[j];
      c22[j] = m1r[j] - m2r[j] + m3r[j] + m6r[j];
    }
  }
}

void StrassenSeq(const double *a_in, std::size_t a_stride_in, const double *b_in, std::size_t b_stride_in, double *c_in,
                 std::size_t c_stride_in, std::size_t n_in) {
  std::function<void(const double *, std::size_t, const double *, std::size_t, double *, std::size_t, std::size_t)>
      impl = [&](const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *c,
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

    std::vector<double> lhs(half * half);
    std::vector<double> rhs(half * half);
    std::vector<double> m1(half * half);
    std::vector<double> m2(half * half);
    std::vector<double> m3(half * half);
    std::vector<double> m4(half * half);
    std::vector<double> m5(half * half);
    std::vector<double> m6(half * half);
    std::vector<double> m7(half * half);

    AddToBuffer(a11, a_stride, a22, a_stride, lhs.data(), half, 1.0);
    AddToBuffer(b11, b_stride, b22, b_stride, rhs.data(), half, 1.0);
    impl(lhs.data(), half, rhs.data(), half, m1.data(), half, half);

    AddToBuffer(a21, a_stride, a22, a_stride, lhs.data(), half, 1.0);
    impl(lhs.data(), half, b11, b_stride, m2.data(), half, half);

    AddToBuffer(b12, b_stride, b22, b_stride, rhs.data(), half, -1.0);
    impl(a11, a_stride, rhs.data(), half, m3.data(), half, half);

    AddToBuffer(b21, b_stride, b11, b_stride, rhs.data(), half, -1.0);
    impl(a22, a_stride, rhs.data(), half, m4.data(), half, half);

    AddToBuffer(a11, a_stride, a12, a_stride, lhs.data(), half, 1.0);
    impl(lhs.data(), half, b22, b_stride, m5.data(), half, half);

    AddToBuffer(a21, a_stride, a11, a_stride, lhs.data(), half, -1.0);
    AddToBuffer(b11, b_stride, b12, b_stride, rhs.data(), half, 1.0);
    impl(lhs.data(), half, rhs.data(), half, m6.data(), half, half);

    AddToBuffer(a12, a_stride, a22, a_stride, lhs.data(), half, -1.0);
    AddToBuffer(b21, b_stride, b22, b_stride, rhs.data(), half, 1.0);
    impl(lhs.data(), half, rhs.data(), half, m7.data(), half, half);

    CombineQuadrants(m1, m2, m3, m4, m5, m6, m7, c, c_stride, half);
  };

  impl(a_in, a_stride_in, b_in, b_stride_in, c_in, c_stride_in, n_in);
}

void StrassenTopStl(const double *a, std::size_t a_stride, const double *b, std::size_t b_stride, double *c,
                    std::size_t c_stride, std::size_t n) {
  if (n <= kCutoff || ppc::util::GetNumThreads() <= 1) {
    StrassenSeq(a, a_stride, b, b_stride, c, c_stride, n);
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
    std::vector<double> lhs(half * half);
    std::vector<double> rhs(half * half);
    AddToBuffer(a11, a_stride, a22, a_stride, lhs.data(), half, 1.0);
    AddToBuffer(b11, b_stride, b22, b_stride, rhs.data(), half, 1.0);
    m1.assign(half * half, 0.0);
    StrassenSeq(lhs.data(), half, rhs.data(), half, m1.data(), half, half);
  });
  auto f2 = std::async(std::launch::async, [&] {
    std::vector<double> lhs(half * half);
    AddToBuffer(a21, a_stride, a22, a_stride, lhs.data(), half, 1.0);
    m2.assign(half * half, 0.0);
    StrassenSeq(lhs.data(), half, b11, b_stride, m2.data(), half, half);
  });
  auto f3 = std::async(std::launch::async, [&] {
    std::vector<double> rhs(half * half);
    AddToBuffer(b12, b_stride, b22, b_stride, rhs.data(), half, -1.0);
    m3.assign(half * half, 0.0);
    StrassenSeq(a11, a_stride, rhs.data(), half, m3.data(), half, half);
  });
  auto f4 = std::async(std::launch::async, [&] {
    std::vector<double> rhs(half * half);
    AddToBuffer(b21, b_stride, b11, b_stride, rhs.data(), half, -1.0);
    m4.assign(half * half, 0.0);
    StrassenSeq(a22, a_stride, rhs.data(), half, m4.data(), half, half);
  });
  auto f5 = std::async(std::launch::async, [&] {
    std::vector<double> lhs(half * half);
    AddToBuffer(a11, a_stride, a12, a_stride, lhs.data(), half, 1.0);
    m5.assign(half * half, 0.0);
    StrassenSeq(lhs.data(), half, b22, b_stride, m5.data(), half, half);
  });
  auto f6 = std::async(std::launch::async, [&] {
    std::vector<double> lhs(half * half);
    std::vector<double> rhs(half * half);
    AddToBuffer(a21, a_stride, a11, a_stride, lhs.data(), half, -1.0);
    AddToBuffer(b11, b_stride, b12, b_stride, rhs.data(), half, 1.0);
    m6.assign(half * half, 0.0);
    StrassenSeq(lhs.data(), half, rhs.data(), half, m6.data(), half, half);
  });
  auto f7 = std::async(std::launch::async, [&] {
    std::vector<double> lhs(half * half);
    std::vector<double> rhs(half * half);
    AddToBuffer(a12, a_stride, a22, a_stride, lhs.data(), half, -1.0);
    AddToBuffer(b21, b_stride, b22, b_stride, rhs.data(), half, 1.0);
    m7.assign(half * half, 0.0);
    StrassenSeq(lhs.data(), half, rhs.data(), half, m7.data(), half, half);
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

}  // namespace

MuhammadkhonIStressenAlgSTL::MuhammadkhonIStressenAlgSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool MuhammadkhonIStressenAlgSTL::ValidationImpl() {
  const auto &in = GetInput();
  return in.a_rows > 0 && in.a_cols_b_rows > 0 && in.b_cols > 0 &&
         in.a.size() == static_cast<size_t>(in.a_rows * in.a_cols_b_rows) &&
         in.b.size() == static_cast<size_t>(in.a_cols_b_rows * in.b_cols);
}

bool MuhammadkhonIStressenAlgSTL::PreProcessingImpl() {
  GetOutput() = {};
  const auto &in = GetInput();
  a_rows_ = in.a_rows;
  a_cols_b_rows_ = in.a_cols_b_rows;
  b_cols_ = in.b_cols;

  size_t max_dim = std::max({a_rows_, a_cols_b_rows_, b_cols_});
  padded_n_ = NextPow2(max_dim);

  padded_a_.assign(padded_n_ * padded_n_, 0.0);
  padded_b_.assign(padded_n_ * padded_n_, 0.0);

  for (size_t i = 0; i < a_rows_; ++i) {
    for (size_t j = 0; j < a_cols_b_rows_; ++j) {
      padded_a_[(i * padded_n_) + j] = in.a[(i * a_cols_b_rows_) + j];
    }
  }

  for (size_t i = 0; i < a_cols_b_rows_; ++i) {
    for (size_t j = 0; j < b_cols_; ++j) {
      padded_b_[(i * padded_n_) + j] = in.b[(i * b_cols_) + j];
    }
  }

  return true;
}

bool MuhammadkhonIStressenAlgSTL::RunImpl() {
  result_c_.assign(padded_n_ * padded_n_, 0.0);

  StrassenTopStl(padded_a_.data(), padded_n_, padded_b_.data(), padded_n_, result_c_.data(), padded_n_, padded_n_);

  auto &out = GetOutput();
  out.assign(a_rows_ * b_cols_, 0.0);
  for (size_t i = 0; i < a_rows_; ++i) {
    for (size_t j = 0; j < b_cols_; ++j) {
      out[(i * b_cols_) + j] = result_c_[(i * padded_n_) + j];
    }
  }

  return true;
}

bool MuhammadkhonIStressenAlgSTL::PostProcessingImpl() {
  return true;
}

}  // namespace muhammadkhon_i_stressen_alg
