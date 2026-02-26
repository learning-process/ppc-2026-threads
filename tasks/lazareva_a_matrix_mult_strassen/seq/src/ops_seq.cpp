#include "lazareva_a_matrix_mult_strassen/seq/include/ops_seq.hpp"

#include <cstddef>
#include <cstdint>
#include <stack>
#include <utility>
#include <vector>

#include "lazareva_a_matrix_mult_strassen/common/include/common.hpp"
#include "task/include/task.hpp"

namespace lazareva_a_matrix_mult_strassen {

namespace {

struct StrassenTask {
  std::vector<double> a;
  std::vector<double> b;
  int n;
  int result_idx;
};

}  // namespace

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
  const auto expected = static_cast<size_t>(n) * static_cast<size_t>(n);
  return std::cmp_equal(GetInput().a.size(), expected) && std::cmp_equal(GetInput().b.size(), expected);
}

bool LazarevaATestTaskSEQ::PreProcessingImpl() {
  n_ = GetInput().n;
  padded_n_ = NextPowerOfTwo(n_);
  a_ = PadMatrix(GetInput().a, n_, padded_n_);
  b_ = PadMatrix(GetInput().b, n_, padded_n_);
  const auto padded_size = static_cast<size_t>(padded_n_) * static_cast<size_t>(padded_n_);
  result_.assign(padded_size, 0.0);
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

std::vector<double> LazarevaATestTaskSEQ::UnpadMatrix(const std::vector<double> &m, int old_n, int new_n) {
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

std::vector<double> LazarevaATestTaskSEQ::Add(const std::vector<double> &a, const std::vector<double> &b, int n) {
  const auto size = static_cast<size_t>(n) * static_cast<size_t>(n);
  std::vector<double> result(size);
  for (size_t i = 0; i < size; ++i) {
    result[i] = a[i] + b[i];
  }
  return result;
}

std::vector<double> LazarevaATestTaskSEQ::Sub(const std::vector<double> &a, const std::vector<double> &b, int n) {
  const auto size = static_cast<size_t>(n) * static_cast<size_t>(n);
  std::vector<double> result(size);
  for (size_t i = 0; i < size; ++i) {
    result[i] = a[i] - b[i];
  }
  return result;
}

void LazarevaATestTaskSEQ::Split(const std::vector<double> &parent, int n, std::vector<double> &a11,
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

std::vector<double> LazarevaATestTaskSEQ::Merge(const std::vector<double> &c11, const std::vector<double> &c12,
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

std::vector<double> LazarevaATestTaskSEQ::NaiveMult(const std::vector<double> &a, const std::vector<double> &b, int n) {
  const auto size = static_cast<size_t>(n) * static_cast<size_t>(n);
  std::vector<double> c(size, 0.0);
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

std::vector<double> LazarevaATestTaskSEQ::Strassen(const std::vector<double> &a, const std::vector<double> &b, int n) {
  struct Frame {
    std::vector<double> a;
    std::vector<double> b;
    int n;
    int phase;
    std::vector<std::vector<double>> m;
    std::vector<double> a11, a12, a21, a22;
    std::vector<double> b11, b12, b21, b22;
  };

  std::stack<Frame> stack;
  std::vector<std::vector<double>> results;

  stack.push({a, b, n, 0, {}, {}, {}, {}, {}, {}, {}, {}, {}});

  while (!stack.empty()) {
    auto &frame = stack.top();

    if (frame.n <= 64) {
      results.push_back(NaiveMult(frame.a, frame.b, frame.n));
      stack.pop();
      continue;
    }

    const int half = frame.n / 2;

    if (frame.phase == 0) {
      Split(frame.a, frame.n, frame.a11, frame.a12, frame.a21, frame.a22);
      Split(frame.b, frame.n, frame.b11, frame.b12, frame.b21, frame.b22);
      frame.m.resize(7);
      frame.phase = 1;

      std::vector<std::pair<std::vector<double>, std::vector<double>>> sub = {
          {Add(frame.a11, frame.a22, half), Add(frame.b11, frame.b22, half)},
          {Add(frame.a21, frame.a22, half), frame.b11},
          {frame.a11, Sub(frame.b12, frame.b22, half)},
          {frame.a22, Sub(frame.b21, frame.b11, half)},
          {Add(frame.a11, frame.a12, half), frame.b22},
          {Sub(frame.a21, frame.a11, half), Add(frame.b11, frame.b12, half)},
          {Sub(frame.a12, frame.a22, half), Add(frame.b21, frame.b22, half)},
      };

      for (int idx = 6; idx >= 0; --idx) {
        stack.push({std::move(sub[static_cast<size_t>(idx)].first),
                    std::move(sub[static_cast<size_t>(idx)].second),
                    half,
                    0,
                    {},
                    {},
                    {},
                    {},
                    {},
                    {},
                    {},
                    {},
                    {}});
      }
      continue;
    }

    const size_t results_needed = 7;
    if (results.size() >= results_needed) {
      const size_t base = results.size() - results_needed;
      for (int idx = 0; idx < 7; ++idx) {
        frame.m[static_cast<size_t>(idx)] = std::move(results[base + static_cast<size_t>(idx)]);
      }
      results.resize(base);

      auto c11 = Add(Sub(Add(frame.m[0], frame.m[3], half), frame.m[4], half), frame.m[6], half);
      auto c12 = Add(frame.m[2], frame.m[4], half);
      auto c21 = Add(frame.m[1], frame.m[3], half);
      auto c22 = Add(Sub(Add(frame.m[0], frame.m[2], half), frame.m[1], half), frame.m[5], half);

      results.push_back(Merge(c11, c12, c21, c22, half));
      stack.pop();
    }
  }

  return results.front();
}

}  // namespace lazareva_a_matrix_mult_strassen
