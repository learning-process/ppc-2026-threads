#include "lazareva_a_matrix_mult_strassen/seq/include/ops_seq.hpp"

#include <vector>

namespace lazareva_a_matrix_mult_strassen {

LazarevaATestTaskSEQ::LazarevaATestTaskSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool LazarevaATestTaskSEQ::ValidationImpl() {
  const int n = GetInput().n;
  if (n <= 0) {
    return false;
  }
  const int expected = n * n;
  return static_cast<int>(GetInput().a.size()) == expected && static_cast<int>(GetInput().b.size()) == expected;
}

bool LazarevaATestTaskSEQ::PreProcessingImpl() {
  n_ = GetInput().n;
  padded_n_ = NextPowerOfTwo(n_);
  a_ = PadMatrix(GetInput().a, n_, padded_n_);
  b_ = PadMatrix(GetInput().b, n_, padded_n_);
  result_.assign(static_cast<size_t>(padded_n_ * padded_n_), 0.0);
  return true;
}

bool LazarevaATestTaskSEQ::RunImpl() {
  result_ = Strassen(a_, b_, padded_n_);
  return true;
}

bool LazarevaATestTaskSEQ::PostProcessingImpl() {
  GetOutput() = UnpadMatrix(result_, padded_n_, n_);
  return true;
}

int LazarevaATestTaskSEQ::NextPowerOfTwo(int n) {
  if (n <= 0) {
    return 1;
  }
  int p = 1;
  while (p < n) {
    p <<= 1;
  }
  return p;
}

std::vector<double> LazarevaATestTaskSEQ::PadMatrix(const std::vector<double> &m, int old_n, int new_n) {
  std::vector<double> result(static_cast<size_t>(new_n * new_n), 0.0);
  for (int i = 0; i < old_n; ++i) {
    for (int j = 0; j < old_n; ++j) {
      result[static_cast<size_t>(i * new_n + j)] = m[static_cast<size_t>(i * old_n + j)];
    }
  }
  return result;
}

std::vector<double> LazarevaATestTaskSEQ::UnpadMatrix(const std::vector<double> &m, int old_n, int new_n) {
  std::vector<double> result(static_cast<size_t>(new_n * new_n));
  for (int i = 0; i < new_n; ++i) {
    for (int j = 0; j < new_n; ++j) {
      result[static_cast<size_t>(i * new_n + j)] = m[static_cast<size_t>(i * old_n + j)];
    }
  }
  return result;
}

std::vector<double> LazarevaATestTaskSEQ::Add(const std::vector<double> &a, const std::vector<double> &b, int n) {
  std::vector<double> result(static_cast<size_t>(n * n));
  for (int i = 0; i < n * n; ++i) {
    result[static_cast<size_t>(i)] = a[static_cast<size_t>(i)] + b[static_cast<size_t>(i)];
  }
  return result;
}

std::vector<double> LazarevaATestTaskSEQ::Sub(const std::vector<double> &a, const std::vector<double> &b, int n) {
  std::vector<double> result(static_cast<size_t>(n * n));
  for (int i = 0; i < n * n; ++i) {
    result[static_cast<size_t>(i)] = a[static_cast<size_t>(i)] - b[static_cast<size_t>(i)];
  }
  return result;
}

void LazarevaATestTaskSEQ::Split(const std::vector<double> &parent, int n, std::vector<double> &a11,
                                 std::vector<double> &a12, std::vector<double> &a21, std::vector<double> &a22) {
  const int half = n / 2;
  a11.resize(static_cast<size_t>(half * half));
  a12.resize(static_cast<size_t>(half * half));
  a21.resize(static_cast<size_t>(half * half));
  a22.resize(static_cast<size_t>(half * half));

  for (int i = 0; i < half; ++i) {
    for (int j = 0; j < half; ++j) {
      const size_t idx = static_cast<size_t>(i * half + j);
      a11[idx] = parent[static_cast<size_t>(i * n + j)];
      a12[idx] = parent[static_cast<size_t>(i * n + j + half)];
      a21[idx] = parent[static_cast<size_t>((i + half) * n + j)];
      a22[idx] = parent[static_cast<size_t>((i + half) * n + j + half)];
    }
  }
}

std::vector<double> LazarevaATestTaskSEQ::Merge(const std::vector<double> &c11, const std::vector<double> &c12,
                                                const std::vector<double> &c21, const std::vector<double> &c22, int n) {
  const int full = n * 2;
  std::vector<double> result(static_cast<size_t>(full * full));

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      const size_t src = static_cast<size_t>(i * n + j);
      result[static_cast<size_t>(i * full + j)] = c11[src];
      result[static_cast<size_t>(i * full + j + n)] = c12[src];
      result[static_cast<size_t>((i + n) * full + j)] = c21[src];
      result[static_cast<size_t>((i + n) * full + j + n)] = c22[src];
    }
  }
  return result;
}

std::vector<double> LazarevaATestTaskSEQ::NaiveMult(const std::vector<double> &a, const std::vector<double> &b, int n) {
  std::vector<double> c(static_cast<size_t>(n * n), 0.0);
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < n; ++k) {
      const double aik = a[static_cast<size_t>(i * n + k)];
      for (int j = 0; j < n; ++j) {
        c[static_cast<size_t>(i * n + j)] += aik * b[static_cast<size_t>(k * n + j)];
      }
    }
  }
  return c;
}

std::vector<double> LazarevaATestTaskSEQ::Strassen(const std::vector<double> &a, const std::vector<double> &b, int n) {
  if (n <= 64) {
    return NaiveMult(a, b, n);
  }

  const int half = n / 2;

  std::vector<double> a11, a12, a21, a22;
  std::vector<double> b11, b12, b21, b22;
  Split(a, n, a11, a12, a21, a22);
  Split(b, n, b11, b12, b21, b22);

  auto m1 = Strassen(Add(a11, a22, half), Add(b11, b22, half), half);
  auto m2 = Strassen(Add(a21, a22, half), b11, half);
  auto m3 = Strassen(a11, Sub(b12, b22, half), half);
  auto m4 = Strassen(a22, Sub(b21, b11, half), half);
  auto m5 = Strassen(Add(a11, a12, half), b22, half);
  auto m6 = Strassen(Sub(a21, a11, half), Add(b11, b12, half), half);
  auto m7 = Strassen(Sub(a12, a22, half), Add(b21, b22, half), half);

  auto c11 = Add(Sub(Add(m1, m4, half), m5, half), m7, half);
  auto c12 = Add(m3, m5, half);
  auto c21 = Add(m2, m4, half);
  auto c22 = Add(Sub(Add(m1, m3, half), m2, half), m6, half);

  return Merge(c11, c12, c21, c22, half);
}

}  // namespace lazareva_a_matrix_mult_strassen
