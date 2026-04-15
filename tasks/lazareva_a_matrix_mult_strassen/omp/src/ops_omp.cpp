#include "lazareva_a_matrix_mult_strassen/omp/include/ops_omp.hpp"

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include "lazareva_a_matrix_mult_strassen/common/include/common.hpp"

namespace lazareva_a_matrix_mult_strassen {

LazarevaATestTaskOMP::LazarevaATestTaskOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool LazarevaATestTaskOMP::ValidationImpl() {
  const int n = GetInput().n;
  if (n <= 0) {
    return false;
  }
  const auto expected = static_cast<size_t>(n) * static_cast<size_t>(n);
  return std::cmp_equal(GetInput().a.size(), expected) && std::cmp_equal(GetInput().b.size(), expected);
}

bool LazarevaATestTaskOMP::PreProcessingImpl() {
  n_ = GetInput().n;
  padded_n_ = NextPowerOfTwo(n_);
  a_ = PadMatrix(GetInput().a, n_, padded_n_);
  b_ = PadMatrix(GetInput().b, n_, padded_n_);
  const auto padded_size = static_cast<size_t>(padded_n_) * static_cast<size_t>(padded_n_);
  result_.assign(padded_size, 0.0);
  return true;
}

bool LazarevaATestTaskOMP::RunImpl() {
  result_ = Strassen(a_, b_, padded_n_);
  return true;
}

bool LazarevaATestTaskOMP::PostProcessingImpl() {
  GetOutput() = UnpadMatrix(result_, padded_n_, n_);
  return true;
}

int LazarevaATestTaskOMP::NextPowerOfTwo(int n) {
  if (n <= 0) {
    return 1;
  }
  int p = 1;
  while (p < n) {
    p <<= 1;
  }
  return p;
}

std::vector<double> LazarevaATestTaskOMP::PadMatrix(const std::vector<double> &m, int old_n, int new_n) {
  const auto new_size = static_cast<size_t>(new_n) * static_cast<size_t>(new_n);
  std::vector<double> result(new_size, 0.0);
  for (int i = 0; i < old_n; ++i) {
    for (int j = 0; j < old_n; ++j) {
      const auto dst = (static_cast<ptrdiff_t>(i) * new_n) + j;
      const auto src = (static_cast<ptrdiff_t>(i) * old_n) + j;
      result[static_cast<size_t>(dst)] = m[static_cast<size_t>(src)];
    }
  }
  return result;
}

std::vector<double> LazarevaATestTaskOMP::UnpadMatrix(const std::vector<double> &m, int old_n, int new_n) {
  const auto new_size = static_cast<size_t>(new_n) * static_cast<size_t>(new_n);
  std::vector<double> result(new_size);
  for (int i = 0; i < new_n; ++i) {
    for (int j = 0; j < new_n; ++j) {
      const auto dst = (static_cast<ptrdiff_t>(i) * new_n) + j;
      const auto src = (static_cast<ptrdiff_t>(i) * old_n) + j;
      result[static_cast<size_t>(dst)] = m[static_cast<size_t>(src)];
    }
  }
  return result;
}

std::vector<double> LazarevaATestTaskOMP::Add(const std::vector<double> &a, const std::vector<double> &b, int n) {
  const auto size = static_cast<size_t>(n) * static_cast<size_t>(n);
  std::vector<double> result(size);
  for (size_t i = 0; i < size; ++i) {
    result[i] = a[i] + b[i];
  }
  return result;
}

std::vector<double> LazarevaATestTaskOMP::Sub(const std::vector<double> &a, const std::vector<double> &b, int n) {
  const auto size = static_cast<size_t>(n) * static_cast<size_t>(n);
  std::vector<double> result(size);
  for (size_t i = 0; i < size; ++i) {
    result[i] = a[i] - b[i];
  }
  return result;
}

void LazarevaATestTaskOMP::Split(const std::vector<double> &parent, int n, std::vector<double> &a11,
                                 std::vector<double> &a12, std::vector<double> &a21, std::vector<double> &a22) {
  const int half = n / 2;
  const auto half_size = static_cast<size_t>(half) * static_cast<size_t>(half);
  a11.resize(half_size);
  a12.resize(half_size);
  a21.resize(half_size);
  a22.resize(half_size);

  for (int i = 0; i < half; ++i) {
    for (int j = 0; j < half; ++j) {
      const auto idx = static_cast<size_t>((static_cast<ptrdiff_t>(i) * half) + j);
      a11[idx] = parent[static_cast<size_t>((static_cast<ptrdiff_t>(i) * n) + j)];
      a12[idx] = parent[static_cast<size_t>((static_cast<ptrdiff_t>(i) * n) + j + half)];
      a21[idx] = parent[static_cast<size_t>((static_cast<ptrdiff_t>(i + half) * n) + j)];
      a22[idx] = parent[static_cast<size_t>((static_cast<ptrdiff_t>(i + half) * n) + j + half)];
    }
  }
}

std::vector<double> LazarevaATestTaskOMP::Merge(const std::vector<double> &c11, const std::vector<double> &c12,
                                                const std::vector<double> &c21, const std::vector<double> &c22, int n) {
  const int full = n * 2;
  const auto full_size = static_cast<size_t>(full) * static_cast<size_t>(full);
  std::vector<double> result(full_size);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      const auto src = static_cast<size_t>((static_cast<ptrdiff_t>(i) * n) + j);
      result[static_cast<size_t>((static_cast<ptrdiff_t>(i) * full) + j)] = c11[src];
      result[static_cast<size_t>((static_cast<ptrdiff_t>(i) * full) + j + n)] = c12[src];
      result[static_cast<size_t>((static_cast<ptrdiff_t>(i + n) * full) + j)] = c21[src];
      result[static_cast<size_t>((static_cast<ptrdiff_t>(i + n) * full) + j + n)] = c22[src];
    }
  }
  return result;
}

std::vector<double> LazarevaATestTaskOMP::NaiveMult(const std::vector<double> &a, const std::vector<double> &b, int n) {
  const auto size = static_cast<size_t>(n) * static_cast<size_t>(n);
  std::vector<double> c(size, 0.0);

#pragma omp parallel for schedule(static) default(none) shared(a, b, c, n)
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < n; ++k) {
      const double aik = a[static_cast<size_t>((static_cast<ptrdiff_t>(i) * n) + k)];
      for (int j = 0; j < n; ++j) {
        c[static_cast<size_t>((static_cast<ptrdiff_t>(i) * n) + j)] +=
            aik * b[static_cast<size_t>((static_cast<ptrdiff_t>(k) * n) + j)];
      }
    }
  }
  return c;
}

std::vector<double> LazarevaATestTaskOMP::Strassen(const std::vector<double> &root_a, const std::vector<double> &root_b,
                                                   int root_n) {
  if (root_n <= 64) {
    return NaiveMult(root_a, root_b, root_n);
  }

  const int half = root_n / 2;

  std::vector<double> a11;
  std::vector<double> a12;
  std::vector<double> a21;
  std::vector<double> a22;
  std::vector<double> b11;
  std::vector<double> b12;
  std::vector<double> b21;
  std::vector<double> b22;
  Split(root_a, root_n, a11, a12, a21, a22);
  Split(root_b, root_n, b11, b12, b21, b22);

  std::array<std::vector<double>, 7> lhs;
  std::array<std::vector<double>, 7> rhs;
  lhs.at(0) = Add(a11, a22, half);
  rhs.at(0) = Add(b11, b22, half);
  lhs.at(1) = Add(a21, a22, half);
  rhs.at(1) = b11;
  lhs.at(2) = a11;
  rhs.at(2) = Sub(b12, b22, half);
  lhs.at(3) = a22;
  rhs.at(3) = Sub(b21, b11, half);
  lhs.at(4) = Add(a11, a12, half);
  rhs.at(4) = b22;
  lhs.at(5) = Sub(a21, a11, half);
  rhs.at(5) = Add(b11, b12, half);
  lhs.at(6) = Sub(a12, a22, half);
  rhs.at(6) = Add(b21, b22, half);

  std::array<std::vector<double>, 7> m;

#pragma omp parallel for schedule(static) default(none) shared(m, lhs, rhs, half)
  for (int k = 0; k < 7; ++k) {
    const auto uk = static_cast<size_t>(k);
    const int nn = half;
    const auto sz = static_cast<size_t>(nn) * static_cast<size_t>(nn);
    std::vector<double> c(sz, 0.0);
    for (int i = 0; i < nn; ++i) {
      for (int ki = 0; ki < nn; ++ki) {
        const double aik = lhs.at(uk)[static_cast<size_t>((static_cast<ptrdiff_t>(i) * nn) + ki)];
        for (int j = 0; j < nn; ++j) {
          c[static_cast<size_t>((static_cast<ptrdiff_t>(i) * nn) + j)] +=
              aik * rhs.at(uk)[static_cast<size_t>((static_cast<ptrdiff_t>(ki) * nn) + j)];
        }
      }
    }
    m.at(uk) = std::move(c);
  }

  auto c11 = Add(Sub(Add(m.at(0), m.at(3), half), m.at(4), half), m.at(6), half);
  auto c12 = Add(m.at(2), m.at(4), half);
  auto c21 = Add(m.at(1), m.at(3), half);
  auto c22 = Add(Sub(Add(m.at(0), m.at(2), half), m.at(1), half), m.at(5), half);

  return Merge(c11, c12, c21, c22, half);
}

}  // namespace lazareva_a_matrix_mult_strassen
