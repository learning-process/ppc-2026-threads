#include "zorin_d_strassen_alg_matrix_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>

namespace zorin_d_strassen_alg_matrix_seq {

namespace {

constexpr std::size_t kCutoff = 64;

inline std::size_t NextPow2(std::size_t x) {
  if (x <= 1) {
    return 1;
  }
  std::size_t p = 1;
  while (p < x) {
    p <<= 1;
  }
  return p;
}

inline void NaiveMul(const std::vector<double> &A, const std::vector<double> &B, std::vector<double> &C,
                     std::size_t n) {
  std::fill(C.begin(), C.end(), 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t k = 0; k < n; ++k) {
      const double aik = A[i * n + k];
      for (std::size_t j = 0; j < n; ++j) {
        C[i * n + j] += aik * B[k * n + j];
      }
    }
  }
}

inline std::vector<double> Add(const std::vector<double> &X, const std::vector<double> &Y) {
  assert(X.size() == Y.size());
  std::vector<double> R(X.size());
  for (std::size_t i = 0; i < X.size(); ++i) {
    R[i] = X[i] + Y[i];
  }
  return R;
}

inline std::vector<double> Sub(const std::vector<double> &X, const std::vector<double> &Y) {
  assert(X.size() == Y.size());
  std::vector<double> R(X.size());
  for (std::size_t i = 0; i < X.size(); ++i) {
    R[i] = X[i] - Y[i];
  }
  return R;
}

inline void Split(const std::vector<double> &M, std::size_t n, std::vector<double> &M11, std::vector<double> &M12,
                  std::vector<double> &M21, std::vector<double> &M22) {
  const std::size_t k = n / 2;
  M11.assign(k * k, 0.0);
  M12.assign(k * k, 0.0);
  M21.assign(k * k, 0.0);
  M22.assign(k * k, 0.0);

  for (std::size_t i = 0; i < k; ++i) {
    for (std::size_t j = 0; j < k; ++j) {
      M11[i * k + j] = M[i * n + j];
      M12[i * k + j] = M[i * n + (j + k)];
      M21[i * k + j] = M[(i + k) * n + j];
      M22[i * k + j] = M[(i + k) * n + (j + k)];
    }
  }
}

inline std::vector<double> Join(const std::vector<double> &C11, const std::vector<double> &C12,
                                const std::vector<double> &C21, const std::vector<double> &C22, std::size_t n) {
  const std::size_t k = n / 2;
  std::vector<double> C(n * n, 0.0);
  for (std::size_t i = 0; i < k; ++i) {
    for (std::size_t j = 0; j < k; ++j) {
      C[i * n + j] = C11[i * k + j];
      C[i * n + (j + k)] = C12[i * k + j];
      C[(i + k) * n + j] = C21[i * k + j];
      C[(i + k) * n + (j + k)] = C22[i * k + j];
    }
  }
  return C;
}

std::vector<double> StrassenRec(const std::vector<double> &A, const std::vector<double> &B, std::size_t n) {
  if (n <= kCutoff) {
    std::vector<double> C(n * n);
    NaiveMul(A, B, C, n);
    return C;
  }

  const std::size_t k = n / 2;
  std::vector<double> A11, A12, A21, A22;
  std::vector<double> B11, B12, B21, B22;
  Split(A, n, A11, A12, A21, A22);
  Split(B, n, B11, B12, B21, B22);

  auto M1 = StrassenRec(Add(A11, A22), Add(B11, B22), k);
  auto M2 = StrassenRec(Add(A21, A22), B11, k);
  auto M3 = StrassenRec(A11, Sub(B12, B22), k);
  auto M4 = StrassenRec(A22, Sub(B21, B11), k);
  auto M5 = StrassenRec(Add(A11, A12), B22, k);
  auto M6 = StrassenRec(Sub(A21, A11), Add(B11, B12), k);
  auto M7 = StrassenRec(Sub(A12, A22), Add(B21, B22), k);

  auto C11 = Add(Sub(Add(M1, M4), M5), M7);
  auto C12 = Add(M3, M5);
  auto C21 = Add(M2, M4);
  auto C22 = Add(Add(Sub(M1, M2), M3), M6);

  return Join(C11, C12, C21, C22, n);
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
  const std::size_t m = NextPow2(n);

  std::vector<double> A_pad(m * m, 0.0);
  std::vector<double> B_pad(m * m, 0.0);

  for (std::size_t i = 0; i < n; ++i) {
    std::copy_n(&in.A[i * n], n, &A_pad[i * m]);
    std::copy_n(&in.B[i * n], n, &B_pad[i * m]);
  }

  auto C_pad = StrassenRec(A_pad, B_pad, m);

  auto &out = GetOutput();
  for (std::size_t i = 0; i < n; ++i) {
    std::copy_n(&C_pad[i * m], n, &out[i * n]);
  }

  return true;
}

bool ZorinDStrassenAlgMatrixSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace zorin_d_strassen_alg_matrix_seq
