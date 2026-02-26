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

  int max_dim = std::max({a_rows_, a_cols_b_rows_, b_cols_});
  padded_n_ = 1;
  while (padded_n_ < max_dim) {
    padded_n_ *= 2;
  }

  padded_a_.assign(padded_n_ * padded_n_, 0.0);
  padded_b_.assign(padded_n_ * padded_n_, 0.0);

  for (int i = 0; i < a_rows_; ++i) {
    for (int j = 0; j < a_cols_b_rows_; ++j) {
      padded_a_[(i * padded_n_) + j] = in.a[(i * a_cols_b_rows_) + j];
    }
  }

  for (int i = 0; i < a_cols_b_rows_; ++i) {
    for (int j = 0; j < b_cols_; ++j) {
      padded_b_[(i * padded_n_) + j] = in.b[(i * b_cols_) + j];
    }
  }
  return true;
}

bool TabalaevAMatrixMulStrassenSEQ::RunImpl() {
  result_c_ = StrassenMultiply(padded_a_, padded_b_, padded_n_);

  auto &out = GetOutput();
  out.assign(a_rows_ * b_cols_, 0.0);

  for (int i = 0; i < a_rows_; ++i) {
    for (int j = 0; j < b_cols_; ++j) {
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

std::vector<double> TabalaevAMatrixMulStrassenSEQ::StrassenMultiply(const std::vector<double> &mat_a,
                                                                    const std::vector<double> &mat_b, int n) {
  std::stack<StrassenFrame> frames;
  std::stack<std::vector<double>> results;

  frames.push({mat_a, mat_b, n, 0});

  while (!frames.empty()) {
    StrassenFrame current = std::move(frames.top());
    frames.pop();

    if (current.stage == 0) {
      if (current.n <= 32) {
        std::vector<double> res(current.n * current.n, 0.0);
        for (int i = 0; i < current.n; ++i) {
          for (int k = 0; k < current.n; ++k) {
            for (int j = 0; j < current.n; ++j) {
              res[(i * current.n) + j] += current.mat_a[(i * current.n) + k] * current.mat_b[(k * current.n) + j];
            }
          }
        }
        results.push(std::move(res));
        continue;
      }

      int h = current.n / 2;
      int sz = h * h;
      std::vector<double> a11(sz);
      std::vector<double> a12(sz);
      std::vector<double> a21(sz);
      std::vector<double> a22(sz);
      std::vector<double> b11(sz);
      std::vector<double> b12(sz);
      std::vector<double> b21(sz);
      std::vector<double> b22(sz);

      for (int i = 0; i < h; ++i) {
        for (int j = 0; j < h; ++j) {
          a11[(i * h) + j] = current.mat_a[(i * current.n) + j];
          a12[(i * h) + j] = current.mat_a[(i * current.n) + j + h];
          a21[(i * h) + j] = current.mat_a[((i + h) * current.n) + j];
          a22[(i * h) + j] = current.mat_a[((i + h) * current.n) + j + h];
          b11[(i * h) + j] = current.mat_b[(i * current.n) + j];
          b12[(i * h) + j] = current.mat_b[(i * current.n) + j + h];
          b21[(i * h) + j] = current.mat_b[((i + h) * current.n) + j];
          b22[(i * h) + j] = current.mat_b[((i + h) * current.n) + j + h];
        }
      }

      frames.push({{}, {}, current.n, 1});

      frames.push({Subtract(a12, a22), Add(b21, b22), h, 0});
      frames.push({Subtract(a21, a11), Add(b11, b12), h, 0});
      frames.push({Add(a11, a12), b22, h, 0});
      frames.push({a22, Subtract(b21, b11), h, 0});
      frames.push({a11, Subtract(b12, b22), h, 0});
      frames.push({Add(a21, a22), b11, h, 0});
      frames.push({Add(a11, a22), Add(b11, b22), h, 0});

    } else {
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

      int h = current.n / 2;
      std::vector<double> res(static_cast<size_t>(current.n) * current.n);
      for (int i = 0; i < h; ++i) {
        for (int j = 0; j < h; ++j) {
          int idx = i * h + j;
          res[(i * current.n) + j] = p1[idx] + p4[idx] - p5[idx] + p7[idx];
          res[(i * current.n) + j + h] = p3[idx] + p5[idx];
          res[((i + h) * current.n) + j] = p2[idx] + p4[idx];
          res[((i + h) * current.n) + j + h] = p1[idx] - p2[idx] + p3[idx] + p6[idx];
        }
      }
      results.push(std::move(res));
    }
  }
  return results.top();
}

}  // namespace tabalaev_a_matrix_mul_strassen
