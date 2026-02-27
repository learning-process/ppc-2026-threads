#pragma once

#include <algorithm>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

namespace korolev_k_matrix_mult::strassen_impl {

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

inline void NaiveMultiply(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c,
                          size_t n) {
  std::ranges::fill(c, 0.0);
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < n; ++k) {
      double a_ik = a[(i * n) + k];
      for (size_t j = 0; j < n; ++j) {
        c[(i * n) + j] += a_ik * b[(k * n) + j];
      }
    }
  }
}

namespace detail {

void AddBlock(const double *a_ptr, const double *b_ptr, size_t ar, size_t ac, size_t br, size_t bc, size_t n, size_t m,
              std::vector<double> &out) {
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < m; ++j) {
      out[(i * m) + j] = a_ptr[((ar + i) * n) + ac + j] + b_ptr[((br + i) * n) + bc + j];
    }
  }
}

void SubBlock(const double *a_ptr, const double *b_ptr, size_t ar, size_t ac, size_t br, size_t bc, size_t n, size_t m,
              std::vector<double> &out) {
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < m; ++j) {
      out[(i * m) + j] = a_ptr[((ar + i) * n) + ac + j] - b_ptr[((br + i) * n) + bc + j];
    }
  }
}

void CopyBlock(const double *a_ptr, size_t ro, size_t co, size_t n, size_t m, std::vector<double> &out) {
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < m; ++j) {
      out[(i * m) + j] = a_ptr[((ro + i) * n) + co + j];
    }
  }
}

void AssembleStrassenResult(std::vector<double> &c, const std::vector<double> &m1, const std::vector<double> &m2,
                            const std::vector<double> &m3, const std::vector<double> &m4, const std::vector<double> &m5,
                            const std::vector<double> &m6, const std::vector<double> &m7, size_t n, size_t m) {
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < m; ++j) {
      c[(i * n) + j] = m1[(i * m) + j] + m4[(i * m) + j] - m5[(i * m) + j] + m7[(i * m) + j];
      c[(i * n) + j + m] = m3[(i * m) + j] + m5[(i * m) + j];
      c[((m + i) * n) + j] = m2[(i * m) + j] + m4[(i * m) + j];
      c[((m + i) * n) + j + m] = m1[(i * m) + j] - m2[(i * m) + j] + m3[(i * m) + j] + m6[(i * m) + j];
    }
  }
}

}  // namespace detail

template <typename ParallelRun>
void StrassenRecurse(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, size_t n,
                     ParallelRun &&parallel_run) {
  if (n <= kStrassenThreshold) {
    NaiveMultiply(a, b, c, n);
    return;
  }

  size_t m = n / 2;
  size_t sz = m * m;

  std::vector<double> m1(sz);
  std::vector<double> m2(sz);
  std::vector<double> m3(sz);
  std::vector<double> m4(sz);
  std::vector<double> m5(sz);
  std::vector<double> m6(sz);
  std::vector<double> m7(sz);

  std::vector<std::function<void()>> tasks = {
      [&]() {
    std::vector<double> t1(sz);
    std::vector<double> t2(sz);
    detail::AddBlock(a.data(), a.data(), 0, 0, m, m, n, m, t1);
    detail::AddBlock(b.data(), b.data(), 0, 0, m, m, n, m, t2);
    StrassenRecurse(t1, t2, m1, m, std::forward<ParallelRun>(parallel_run));
  },
      [&]() {
    std::vector<double> t1(sz);
    std::vector<double> t2(sz);
    detail::AddBlock(a.data(), a.data(), m, 0, m, m, n, m, t1);
    detail::CopyBlock(b.data(), 0, 0, n, m, t2);
    StrassenRecurse(t1, t2, m2, m, std::forward<ParallelRun>(parallel_run));
  },
      [&]() {
    std::vector<double> t1(sz);
    std::vector<double> t2(sz);
    detail::SubBlock(b.data(), b.data(), 0, m, m, m, n, m, t1);
    detail::CopyBlock(a.data(), 0, 0, n, m, t2);
    StrassenRecurse(t2, t1, m3, m, std::forward<ParallelRun>(parallel_run));
  },
      [&]() {
    std::vector<double> t1(sz);
    std::vector<double> t2(sz);
    detail::SubBlock(b.data(), b.data(), m, 0, 0, 0, n, m, t1);
    detail::CopyBlock(a.data(), m, m, n, m, t2);
    StrassenRecurse(t2, t1, m4, m, std::forward<ParallelRun>(parallel_run));
  },
      [&]() {
    std::vector<double> t1(sz);
    std::vector<double> t2(sz);
    detail::AddBlock(a.data(), a.data(), 0, 0, 0, m, n, m, t1);
    detail::CopyBlock(b.data(), m, m, n, m, t2);
    StrassenRecurse(t1, t2, m5, m, std::forward<ParallelRun>(parallel_run));
  },
      [&]() {
    std::vector<double> t1(sz);
    std::vector<double> t2(sz);
    detail::SubBlock(a.data(), a.data(), m, 0, 0, 0, n, m, t1);
    detail::AddBlock(b.data(), b.data(), 0, 0, 0, m, n, m, t2);
    StrassenRecurse(t1, t2, m6, m, std::forward<ParallelRun>(parallel_run));
  },
      [&]() {
    std::vector<double> t1(sz);
    std::vector<double> t2(sz);
    detail::SubBlock(a.data(), a.data(), 0, m, m, m, n, m, t1);
    detail::AddBlock(b.data(), b.data(), m, 0, m, m, n, m, t2);
    StrassenRecurse(t1, t2, m7, m, std::forward<ParallelRun>(parallel_run));
  },
  };

  parallel_run(tasks);

  detail::AssembleStrassenResult(c, m1, m2, m3, m4, m5, m6, m7, n, m);
}

inline void StrassenMultiply(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c,
                             size_t n, const std::function<void(std::vector<std::function<void()>> &)> &parallel_run) {
  StrassenRecurse(a, b, c, n, parallel_run);
}

}  // namespace korolev_k_matrix_mult::strassen_impl
