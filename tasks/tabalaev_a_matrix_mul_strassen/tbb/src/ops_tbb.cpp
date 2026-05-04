#include "tabalaev_a_matrix_mul_strassen/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "tabalaev_a_matrix_mul_strassen/common/include/common.hpp"

namespace tabalaev_a_matrix_mul_strassen {

static constexpr size_t kBaseCaseSize = 128;
static constexpr size_t kParallelThreshold = 262144;

TabalaevAMatrixMulStrassenTBB::TabalaevAMatrixMulStrassenTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool TabalaevAMatrixMulStrassenTBB::ValidationImpl() {
  const auto &in = GetInput();
  return (in.a_rows > 0) && (in.a_cols_b_rows > 0) && (in.b_cols > 0) &&
         (in.a.size() == static_cast<size_t>(in.a_rows * in.a_cols_b_rows)) &&
         (in.b.size() == static_cast<size_t>(in.a_cols_b_rows * in.b_cols));
}

bool TabalaevAMatrixMulStrassenTBB::PreProcessingImpl() {
  const auto &in = GetInput();

  a_rows_ = in.a_rows;
  a_cols_b_rows_ = in.a_cols_b_rows;
  b_cols_ = in.b_cols;

  size_t max_dim = std::max({a_rows_, a_cols_b_rows_, b_cols_});

  padded_n_ = 1;
  while (padded_n_ < max_dim) {
    padded_n_ *= 2;
  }

  padded_a_.assign(padded_n_ * padded_n_, 0.0);
  padded_b_.assign(padded_n_ * padded_n_, 0.0);

  const size_t n = padded_n_;

  tbb::parallel_for(static_cast<size_t>(0), a_rows_, [this, n, &in](size_t i) {
    for (size_t j = 0; j < a_cols_b_rows_; ++j) {
      padded_a_[(i * n) + j] = in.a[(i * a_cols_b_rows_) + j];
    }
  });

  tbb::parallel_for(static_cast<size_t>(0), a_cols_b_rows_, [this, n, &in](size_t i) {
    for (size_t j = 0; j < b_cols_; ++j) {
      padded_b_[(i * n) + j] = in.b[(i * b_cols_) + j];
    }
  });

  return true;
}

bool TabalaevAMatrixMulStrassenTBB::RunImpl() {
  result_c_ = StrassenMultiply(padded_a_, padded_b_, padded_n_);

  auto &out = GetOutput();
  out.assign(a_rows_ * b_cols_, 0.0);

  const size_t n = padded_n_;

  tbb::parallel_for(static_cast<size_t>(0), a_rows_, [this, n, &out](size_t i) {
    for (size_t j = 0; j < b_cols_; ++j) {
      out[(i * b_cols_) + j] = result_c_[(i * n) + j];
    }
  });

  return true;
}

bool TabalaevAMatrixMulStrassenTBB::PostProcessingImpl() {
  return true;
}

void TabalaevAMatrixMulStrassenTBB::Add(const std::vector<double> &mat_a, const std::vector<double> &mat_b,
                                        std::vector<double> &res) {
  const size_t n = mat_a.size();
  res.resize(n);

  if (n >= kParallelThreshold) {
    tbb::parallel_for(static_cast<size_t>(0), n, [&](size_t i) { res[i] = mat_a[i] + mat_b[i]; });
  } else {
    for (size_t i = 0; i < n; ++i) {
      res[i] = mat_a[i] + mat_b[i];
    }
  }
}

void TabalaevAMatrixMulStrassenTBB::Subtract(const std::vector<double> &mat_a, const std::vector<double> &mat_b,
                                             std::vector<double> &res) {
  const size_t n = mat_a.size();
  res.resize(n);

  if (n >= kParallelThreshold) {
    tbb::parallel_for(static_cast<size_t>(0), n, [&](size_t i) { res[i] = mat_a[i] - mat_b[i]; });
  } else {
    for (size_t i = 0; i < n; ++i) {
      res[i] = mat_a[i] - mat_b[i];
    }
  }
}

std::vector<double> TabalaevAMatrixMulStrassenTBB::BaseMultiply(const std::vector<double> &mat_a,
                                                                const std::vector<double> &mat_b, size_t n) {
  std::vector<double> res(n * n, 0.0);

  tbb::parallel_for(static_cast<size_t>(0), n, [&](size_t i) {
    for (size_t k = 0; k < n; ++k) {
      double temp = mat_a[(i * n) + k];
      if (temp == 0.0) {
        continue;
      }
      for (size_t j = 0; j < n; ++j) {
        res[(i * n) + j] += temp * mat_b[(k * n) + j];
      }
    }
  });

  return res;
}

void TabalaevAMatrixMulStrassenTBB::SplitMatrix(const std::vector<double> &src, size_t n, std::vector<double> &c11,
                                                std::vector<double> &c12, std::vector<double> &c21,
                                                std::vector<double> &c22) {
  size_t h = n / 2;
  size_t sz = h * h;

  c11.resize(sz);
  c12.resize(sz);
  c21.resize(sz);
  c22.resize(sz);

  if (n * n >= kParallelThreshold) {
    tbb::parallel_for(static_cast<size_t>(0), h, [&](size_t i) {
      for (size_t j = 0; j < h; ++j) {
        size_t src_idx = (i * n) + j;
        size_t dst_idx = (i * h) + j;

        c11[dst_idx] = src[src_idx];
        c12[dst_idx] = src[src_idx + h];
        c21[dst_idx] = src[src_idx + (h * n)];
        c22[dst_idx] = src[src_idx + (h * n) + h];
      }
    });
  } else {
    for (size_t i = 0; i < h; ++i) {
      for (size_t j = 0; j < h; ++j) {
        size_t src_idx = (i * n) + j;
        size_t dst_idx = (i * h) + j;

        c11[dst_idx] = src[src_idx];
        c12[dst_idx] = src[src_idx + h];
        c21[dst_idx] = src[src_idx + (h * n)];
        c22[dst_idx] = src[src_idx + (h * n) + h];
      }
    }
  }
}

std::vector<double> TabalaevAMatrixMulStrassenTBB::CombineMatrix(const std::vector<double> &c11,
                                                                 const std::vector<double> &c12,
                                                                 const std::vector<double> &c21,
                                                                 const std::vector<double> &c22, size_t n) {
  size_t h = n / 2;

  std::vector<double> res(n * n);

  if (n * n >= kParallelThreshold) {
    tbb::parallel_for(static_cast<size_t>(0), h, [&](size_t i) {
      for (size_t j = 0; j < h; ++j) {
        size_t src_idx = (i * h) + j;

        res[(i * n) + j] = c11[src_idx];
        res[(i * n) + j + h] = c12[src_idx];
        res[((i + h) * n) + j] = c21[src_idx];
        res[((i + h) * n) + j + h] = c22[src_idx];
      }
    });
  } else {
    for (size_t i = 0; i < h; ++i) {
      for (size_t j = 0; j < h; ++j) {
        size_t src_idx = (i * h) + j;

        res[(i * n) + j] = c11[src_idx];
        res[(i * n) + j + h] = c12[src_idx];
        res[((i + h) * n) + j] = c21[src_idx];
        res[((i + h) * n) + j + h] = c22[src_idx];
      }
    }
  }

  return res;
}

std::vector<double> TabalaevAMatrixMulStrassenTBB::StrassenMultiply(const std::vector<double> &mat_a,
                                                                    const std::vector<double> &mat_b, size_t n) {
  std::vector<StrassenFrameTBB> frames;
  std::vector<std::vector<double>> results;

  frames.reserve(64);
  results.reserve(64);

  frames.push_back({mat_a, mat_b, n, 0});

  std::vector<double> t1;
  std::vector<double> t2;

  while (!frames.empty()) {
    StrassenFrameTBB current = std::move(frames.back());
    frames.pop_back();

    if (current.n <= kBaseCaseSize) {
      results.push_back(BaseMultiply(current.mat_a, current.mat_b, current.n));
      continue;
    }

    if (current.stage == 8) {
      std::vector<std::vector<double>> p(7);

      for (int i = 6; i >= 0; --i) {
        p[i] = std::move(results.back());
        results.pop_back();
      }

      size_t h = current.n / 2;
      size_t sz = h * h;

      std::vector<double> c11(sz);
      std::vector<double> c12(sz);
      std::vector<double> c21(sz);
      std::vector<double> c22(sz);

      tbb::parallel_for(static_cast<size_t>(0), sz, [&](size_t i) {
        c11[i] = p[0][i] + p[3][i] - p[4][i] + p[6][i];
        c12[i] = p[2][i] + p[4][i];
        c21[i] = p[1][i] + p[3][i];
        c22[i] = p[0][i] - p[1][i] + p[2][i] + p[5][i];
      });

      results.push_back(CombineMatrix(c11, c12, c21, c22, current.n));
    } else {
      size_t h = current.n / 2;

      std::vector<double> a11;
      std::vector<double> a12;
      std::vector<double> a21;
      std::vector<double> a22;
      std::vector<double> b11;
      std::vector<double> b12;
      std::vector<double> b21;
      std::vector<double> b22;

      SplitMatrix(current.mat_a, current.n, a11, a12, a21, a22);
      SplitMatrix(current.mat_b, current.n, b11, b12, b21, b22);

      frames.push_back({{}, {}, current.n, 8});

      Subtract(a12, a22, t1);
      Add(b21, b22, t2);
      frames.push_back({t1, t2, h, 0});

      Subtract(a21, a11, t1);
      Add(b11, b12, t2);
      frames.push_back({t1, t2, h, 0});

      Add(a11, a12, t1);
      frames.push_back({t1, b22, h, 0});

      Subtract(b21, b11, t2);
      frames.push_back({a22, t2, h, 0});

      Subtract(b12, b22, t2);
      frames.push_back({a11, t2, h, 0});

      Add(a21, a22, t1);
      frames.push_back({t1, b11, h, 0});

      Add(a11, a22, t1);
      Add(b11, b22, t2);
      frames.push_back({t1, t2, h, 0});
    }
  }

  return std::move(results.back());
}

}  // namespace tabalaev_a_matrix_mul_strassen
