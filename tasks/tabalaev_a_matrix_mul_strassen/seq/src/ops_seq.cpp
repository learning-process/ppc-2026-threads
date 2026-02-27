#include "tabalaev_a_matrix_mul_strassen/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stack>
#include <utility>
#include <vector>

#include "tabalaev_a_matrix_mul_strassen/common/include/common.hpp"

namespace tabalaev_a_matrix_mul_strassen {

TabalaevAMatrixMulStrassenSEQ::TabalaevAMatrixMulStrassenSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool TabalaevAMatrixMulStrassenSEQ::ValidationImpl() {
  const auto &in = GetInput();
  return in.a_rows > 0 && in.a_cols_b_rows > 0 && in.b_cols > 0 &&
         in.a.size() == static_cast<size_t>(in.a_rows * in.a_cols_b_rows) &&
         in.b.size() == static_cast<size_t>(in.a_cols_b_rows * in.b_cols);
}

bool TabalaevAMatrixMulStrassenSEQ::PreProcessingImpl() {
  GetOutput() = {};
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

bool TabalaevAMatrixMulStrassenSEQ::RunImpl() {
  result_c_ = StrassenMultiply(padded_a_, padded_b_, padded_n_);

  auto &out = GetOutput();
  out.assign(a_rows_ * b_cols_, 0.0);

  for (size_t i = 0; i < a_rows_; ++i) {
    for (size_t j = 0; j < b_cols_; ++j) {
      out[(i * b_cols_) + j] = result_c_[(i * padded_n_) + j];
    }
  }
  return true;
}

bool TabalaevAMatrixMulStrassenSEQ::PostProcessingImpl() {
  return true;
}

std::vector<double> TabalaevAMatrixMulStrassenSEQ::Add(const std::vector<double> &mat_a,
                                                       const std::vector<double> &mat_b) {
  std::vector<double> res(mat_a.size());
  for (size_t i = 0; i < mat_a.size(); ++i) {
    res[i] = mat_a[i] + mat_b[i];
  }
  return res;
}

std::vector<double> TabalaevAMatrixMulStrassenSEQ::Subtract(const std::vector<double> &mat_a,
                                                            const std::vector<double> &mat_b) {
  std::vector<double> res(mat_a.size());
  for (size_t i = 0; i < mat_a.size(); ++i) {
    res[i] = mat_a[i] - mat_b[i];
  }
  return res;
}

std::vector<double> TabalaevAMatrixMulStrassenSEQ::BaseMultiply(const std::vector<double> &mat_a,
                                                                const std::vector<double> &mat_b, size_t n) {
  std::vector<double> res(n * n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < n; ++k) {
      for (size_t j = 0; j < n; ++j) {
        res[(i * n) + j] += mat_a[(i * n) + k] * mat_b[(k * n) + j];
      }
    }
  }
  return res;
}

void TabalaevAMatrixMulStrassenSEQ::PushStrassenSubtasks(std::stack<StrassenFrame> &frames,
                                                         const std::vector<double> &mat_a,
                                                         const std::vector<double> &mat_b, size_t n) {
  size_t h = n / 2;
  size_t sz = h * h;
  std::vector<double> a11(sz);
  std::vector<double> a12(sz);
  std::vector<double> a21(sz);
  std::vector<double> a22(sz);
  std::vector<double> b11(sz);
  std::vector<double> b12(sz);
  std::vector<double> b21(sz);
  std::vector<double> b22(sz);

  for (size_t i = 0; i < h; ++i) {
    for (size_t j = 0; j < h; ++j) {
      size_t idx_src = (i * n) + j;
      size_t idx_dst = (i * h) + j;
      a11[idx_dst] = mat_a[idx_src];
      a12[idx_dst] = mat_a[idx_src + h];
      a21[idx_dst] = mat_a[idx_src + (h * n)];
      a22[idx_dst] = mat_a[idx_src + (h * n) + h];

      b11[idx_dst] = mat_b[idx_src];
      b12[idx_dst] = mat_b[idx_src + h];
      b21[idx_dst] = mat_b[idx_src + (h * n)];
      b22[idx_dst] = mat_b[idx_src + (h * n) + h];
    }
  }

  frames.push({{}, {}, n, 1});

  frames.push({Subtract(a12, a22), Add(b21, b22), h, 0});
  frames.push({Subtract(a21, a11), Add(b11, b12), h, 0});
  frames.push({Add(a11, a12), b22, h, 0});
  frames.push({a22, Subtract(b21, b11), h, 0});
  frames.push({a11, Subtract(b12, b22), h, 0});
  frames.push({Add(a21, a22), b11, h, 0});
  frames.push({Add(a11, a22), Add(b11, b22), h, 0});
}

std::vector<double> TabalaevAMatrixMulStrassenSEQ::CombineStrassenResults(std::stack<std::vector<double>> &results,
                                                                          size_t n) {
  auto p7 = std::move(results.top());
  results.pop();
  auto p6 = std::move(results.top());
  results.pop();
  auto p5 = std::move(results.top());
  results.pop();
  auto p4 = std::move(results.top());
  results.pop();
  auto p3 = std::move(results.top());
  results.pop();
  auto p2 = std::move(results.top());
  results.pop();
  auto p1 = std::move(results.top());
  results.pop();

  size_t h = n / 2;
  std::vector<double> res(n * n);
  for (int i = 0; i < h; ++i) {
    for (int j = 0; j < h; ++j) {
      int idx = (i * h) + j;
      res[(i * n) + j] = p1[idx] + p4[idx] - p5[idx] + p7[idx];
      res[(i * n) + j + h] = p3[idx] + p5[idx];
      res[((i + h) * n) + j] = p2[idx] + p4[idx];
      res[((i + h) * n) + j + h] = p1[idx] - p2[idx] + p3[idx] + p6[idx];
    }
  }
  return res;
}

std::vector<double> TabalaevAMatrixMulStrassenSEQ::StrassenMultiply(const std::vector<double> &mat_a,
                                                                    const std::vector<double> &mat_b, size_t n) {
  std::stack<StrassenFrame> frames;
  std::stack<std::vector<double>> results;

  frames.push({mat_a, mat_b, n, 0});

  while (!frames.empty()) {
    StrassenFrame current = std::move(frames.top());
    frames.pop();

    if (current.stage == 0) {
      if (current.n <= 32) {
        results.push(BaseMultiply(current.mat_a, current.mat_b, current.n));
      } else {
        PushStrassenSubtasks(frames, current.mat_a, current.mat_b, current.n);
      }
    } else {
      results.push(CombineStrassenResults(results, current.n));
    }
  }
  return results.top();
}
}  // namespace tabalaev_a_matrix_mul_strassen
