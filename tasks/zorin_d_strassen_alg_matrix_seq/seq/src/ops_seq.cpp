#include "zorin_d_strassen_alg_matrix_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <ranges>
#include <utility>
#include <vector>

namespace zorin_d_strassen_alg_matrix_seq {

namespace {

constexpr std::size_t kCutoff = 64;

std::size_t NextPow2(std::size_t x) {
  if (x <= 1) {
    return 1;
  }
  std::size_t p = 1;
  while (p < x) {
    p <<= 1U;
  }
  return p;
}

void NaiveMul(const std::vector<double> &A, const std::vector<double> &B, std::vector<double> &C, std::size_t n) {
  std::ranges::fill(C, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t iRow = i * n;
    for (std::size_t k = 0; k < n; ++k) {
      const double aik = A[iRow + k];
      const std::size_t kRow = k * n;
      for (std::size_t j = 0; j < n; ++j) {
        C[iRow + j] += aik * B[kRow + j];
      }
    }
  }
}

std::vector<double> AddVec(const std::vector<double> &X, const std::vector<double> &Y) {
  assert(X.size() == Y.size());
  std::vector<double> R(X.size());
  for (std::size_t i = 0; i < X.size(); ++i) {
    R[i] = X[i] + Y[i];
  }
  return R;
}

std::vector<double> SubVec(const std::vector<double> &X, const std::vector<double> &Y) {
  assert(X.size() == Y.size());
  std::vector<double> R(X.size());
  for (std::size_t i = 0; i < X.size(); ++i) {
    R[i] = X[i] - Y[i];
  }
  return R;
}

void SplitMat(const std::vector<double> &M, std::size_t n, std::vector<double> &M11, std::vector<double> &M12,
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

std::vector<double> JoinMat(const std::vector<double> &C11, const std::vector<double> &C12,
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

struct Frame {
  std::size_t n{};
  int parentSlot{-1};

  std::vector<double> A;
  std::vector<double> B;
  std::vector<double> C;

  std::size_t k{};
  int step{0};

  std::vector<double> A11, A12, A21, A22;
  std::vector<double> B11, B12, B21, B22;

  std::vector<double> M1, M2, M3, M4, M5, M6, M7;

  Frame(std::vector<double> a, std::vector<double> b, std::size_t n_, int parentSlot_)
      : n(n_), parentSlot(parentSlot_), A(std::move(a)), B(std::move(b)), C(n_ * n_, 0.0) {}
};

void StoreM(Frame &f, int slot, std::vector<double> v) {
  switch (slot) {
    case 1:
      f.M1 = std::move(v);
      break;
    case 2:
      f.M2 = std::move(v);
      break;
    case 3:
      f.M3 = std::move(v);
      break;
    case 4:
      f.M4 = std::move(v);
      break;
    case 5:
      f.M5 = std::move(v);
      break;
    case 6:
      f.M6 = std::move(v);
      break;
    case 7:
      f.M7 = std::move(v);
      break;
    default:
      break;
  }
}

Frame MakeChildForStep(const Frame &f, int step) {
  const std::size_t k = f.k;

  if (step == 0) {
    return Frame(AddVec(f.A11, f.A22), AddVec(f.B11, f.B22), k, 1);
  }
  if (step == 1) {
    return Frame(AddVec(f.A21, f.A22), f.B11, k, 2);
  }
  if (step == 2) {
    return Frame(f.A11, SubVec(f.B12, f.B22), k, 3);
  }
  if (step == 3) {
    return Frame(f.A22, SubVec(f.B21, f.B11), k, 4);
  }
  if (step == 4) {
    return Frame(AddVec(f.A11, f.A12), f.B22, k, 5);
  }
  if (step == 5) {
    return Frame(SubVec(f.A21, f.A11), AddVec(f.B11, f.B12), k, 6);
  }
  return Frame(SubVec(f.A12, f.A22), AddVec(f.B21, f.B22), k, 7);
}

std::vector<double> StrassenIter(std::vector<double> A, std::vector<double> B, std::size_t n) {
  std::vector<Frame> st;
  st.emplace_back(std::move(A), std::move(B), n, -1);

  std::vector<double> rootOut;

  while (!st.empty()) {
    Frame &f = st.back();

    if (f.n <= kCutoff) {
      NaiveMul(f.A, f.B, f.C, f.n);

      Frame finished = std::move(f);
      st.pop_back();

      if (st.empty()) {
        rootOut = std::move(finished.C);
      } else {
        Frame &parent = st.back();
        StoreM(parent, finished.parentSlot, std::move(finished.C));
      }
      continue;
    }

    if (f.step == 0) {
      f.k = f.n / 2;
      SplitMat(f.A, f.n, f.A11, f.A12, f.A21, f.A22);
      SplitMat(f.B, f.n, f.B11, f.B12, f.B21, f.B22);
    }

    if (f.step < 7) {
      Frame child = MakeChildForStep(f, f.step);
      ++f.step;
      st.emplace_back(std::move(child));
      continue;
    }

    const auto C11 = AddVec(SubVec(AddVec(f.M1, f.M4), f.M5), f.M7);
    const auto C12 = AddVec(f.M3, f.M5);
    const auto C21 = AddVec(f.M2, f.M4);
    const auto C22 = AddVec(AddVec(SubVec(f.M1, f.M2), f.M3), f.M6);

    f.C = JoinMat(C11, C12, C21, C22, f.n);

    Frame finished = std::move(f);
    st.pop_back();

    if (st.empty()) {
      rootOut = std::move(finished.C);
    } else {
      Frame &parent = st.back();
      StoreM(parent, finished.parentSlot, std::move(finished.C));
    }
  }

  return rootOut;
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
  const std::size_t need = in.n * in.n;
  if (in.A.size() != need) {
    return false;
  }
  if (in.B.size() != need) {
    return false;
  }
  if (!GetOutput().empty()) {
    return false;
  }
  return true;
}

bool ZorinDStrassenAlgMatrixSEQ::PreProcessingImpl() {
  const std::size_t n = GetInput().n;
  GetOutput().assign(n * n, 0.0);
  return true;
}

bool ZorinDStrassenAlgMatrixSEQ::RunImpl() {
  const auto &in = GetInput();
  const std::size_t n = in.n;
  const std::size_t m = NextPow2(n);

  std::vector<double> APad(m * m, 0.0);
  std::vector<double> BPad(m * m, 0.0);

  for (std::size_t i = 0; i < n; ++i) {
    std::copy_n(&in.A[i * n], n, &APad[i * m]);
    std::copy_n(&in.B[i * n], n, &BPad[i * m]);
  }

  const auto CPad = StrassenIter(std::move(APad), std::move(BPad), m);

  auto &out = GetOutput();
  for (std::size_t i = 0; i < n; ++i) {
    std::copy_n(&CPad[i * m], n, &out[i * n]);
  }

  return true;
}

bool ZorinDStrassenAlgMatrixSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace zorin_d_strassen_alg_matrix_seq
