#pragma once

#include <algorithm>
#include <cstddef>
#include <functional>
#include <vector>

namespace korolev_k_matrix_mult {
namespace strassen_impl {

constexpr size_t kStrassenThreshold = 64;

inline size_t NextPowerOf2(size_t n) {
  if (n <= 1) {
    return 1;
  }
  size_t p = 1;
  while (p < n) {
    p *= 2;
  }
  return p;
}

inline void NaiveMultiply(const std::vector<double> &A, const std::vector<double> &B, std::vector<double> &C,
                          size_t n) {
  std::fill(C.begin(), C.end(), 0.0);
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < n; ++k) {
      double a_ik = A[i * n + k];
      for (size_t j = 0; j < n; ++j) {
        C[i * n + j] += a_ik * B[k * n + j];
      }
    }
  }
}

template <typename ParallelRun>
void StrassenRecurse(const std::vector<double> &A, const std::vector<double> &B, std::vector<double> &C, size_t n,
                     ParallelRun &&parallel_run) {
  if (n <= kStrassenThreshold) {
    NaiveMultiply(A, B, C, n);
    return;
  }

  size_t m = n / 2;
  size_t sz = m * m;

  std::vector<double> M1(sz), M2(sz), M3(sz), M4(sz), M5(sz), M6(sz), M7(sz);

  auto add_block = [&](const double *a, const double *b, size_t ar, size_t ac, size_t br, size_t bc,
                       std::vector<double> &out) {
    for (size_t i = 0; i < m; ++i) {
      for (size_t j = 0; j < m; ++j) {
        out[i * m + j] = a[(ar + i) * n + ac + j] + b[(br + i) * n + bc + j];
      }
    }
  };
  auto sub_block = [&](const double *a, const double *b, size_t ar, size_t ac, size_t br, size_t bc,
                       std::vector<double> &out) {
    for (size_t i = 0; i < m; ++i) {
      for (size_t j = 0; j < m; ++j) {
        out[i * m + j] = a[(ar + i) * n + ac + j] - b[(br + i) * n + bc + j];
      }
    }
  };
  auto copy_block = [&](const double *a, size_t ro, size_t co, std::vector<double> &out) {
    for (size_t i = 0; i < m; ++i) {
      for (size_t j = 0; j < m; ++j) {
        out[i * m + j] = a[(ro + i) * n + co + j];
      }
    }
  };

  std::vector<std::function<void()>> tasks = {
      [&]() {
    std::vector<double> t1(sz), t2(sz);
    add_block(A.data(), A.data(), 0, 0, m, m, t1);
    add_block(B.data(), B.data(), 0, 0, m, m, t2);
    StrassenRecurse(t1, t2, M1, m, parallel_run);
  },
      [&]() {
    std::vector<double> t1(sz), t2(sz);
    add_block(A.data(), A.data(), m, 0, m, m, t1);
    copy_block(B.data(), 0, 0, t2);
    StrassenRecurse(t1, t2, M2, m, parallel_run);
  },
      [&]() {
    std::vector<double> t1(sz), t2(sz);
    sub_block(B.data(), B.data(), 0, m, m, m, t1);
    copy_block(A.data(), 0, 0, t2);
    StrassenRecurse(t2, t1, M3, m, parallel_run);
  },
      [&]() {
    std::vector<double> t1(sz), t2(sz);
    sub_block(B.data(), B.data(), m, 0, 0, 0, t1);
    copy_block(A.data(), m, m, t2);
    StrassenRecurse(t2, t1, M4, m, parallel_run);
  },
      [&]() {
    std::vector<double> t1(sz), t2(sz);
    add_block(A.data(), A.data(), 0, 0, 0, m, t1);
    copy_block(B.data(), m, m, t2);
    StrassenRecurse(t1, t2, M5, m, parallel_run);
  },
      [&]() {
    std::vector<double> t1(sz), t2(sz);
    sub_block(A.data(), A.data(), m, 0, 0, 0, t1);
    add_block(B.data(), B.data(), 0, 0, 0, m, t2);
    StrassenRecurse(t1, t2, M6, m, parallel_run);
  },
      [&]() {
    std::vector<double> t1(sz), t2(sz);
    sub_block(A.data(), A.data(), 0, m, m, m, t1);
    add_block(B.data(), B.data(), m, 0, m, m, t2);
    StrassenRecurse(t1, t2, M7, m, parallel_run);
  },
  };

  parallel_run(tasks);

  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < m; ++j) {
      C[i * n + j] = M1[i * m + j] + M4[i * m + j] - M5[i * m + j] + M7[i * m + j];
      C[i * n + j + m] = M3[i * m + j] + M5[i * m + j];
      C[(m + i) * n + j] = M2[i * m + j] + M4[i * m + j];
      C[(m + i) * n + j + m] = M1[i * m + j] - M2[i * m + j] + M3[i * m + j] + M6[i * m + j];
    }
  }
}

inline void StrassenMultiply(const std::vector<double> &A, const std::vector<double> &B, std::vector<double> &C,
                             size_t n, const std::function<void(std::vector<std::function<void()>> &)> &parallel_run) {
  StrassenRecurse(A, B, C, n, parallel_run);
}

}  // namespace strassen_impl
}  // namespace korolev_k_matrix_mult
