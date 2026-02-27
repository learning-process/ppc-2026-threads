#include "muhammadkhon_i_stressen_alg/seq/include/ops_seq.hpp"

#include <array>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "muhammadkhon_i_stressen_alg/common/include/common.hpp"

namespace muhammadkhon_i_stressen_alg {

MuhammadkhonIStressenAlgSEQ::MuhammadkhonIStressenAlgSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool MuhammadkhonIStressenAlgSEQ::ValidationImpl() {
  const int n = GetInput().n;
  if (n <= 0) {
    return false;
  }
  const auto expected = static_cast<size_t>(n) * static_cast<size_t>(n);

  if (expected == 0 || expected > GetInput().a.max_size()) {
    return false;
  }

  return std::cmp_equal(GetInput().a.size(), expected) && std::cmp_equal(GetInput().b.size(), expected);
}

bool MuhammadkhonIStressenAlgSEQ::PreProcessingImpl() {
  n_ = GetInput().n;
  padded_n_ = NextPowerOfTwo(n_);
  a_ = PadMatrix(GetInput().a, n_, padded_n_);
  b_ = PadMatrix(GetInput().b, n_, padded_n_);
  const auto padded_size = static_cast<size_t>(padded_n_) * static_cast<size_t>(padded_n_);
  result_.assign(padded_size, 0.0);
  return true;
}

bool MuhammadkhonIStressenAlgSEQ::RunImpl() {
  result_ = Strassen(a_, b_, padded_n_);
  return true;
}

bool MuhammadkhonIStressenAlgSEQ::PostProcessingImpl() {
  GetOutput() = UnpadMatrix(result_, padded_n_, n_);
  return true;
}

std::vector<double> MuhammadkhonIStressenAlgSEQ::PadMatrix(const std::vector<double> &m, int old_n, int new_n) {
  if (old_n > new_n) {
    throw std::runtime_error("new_n must be >= old_n for padding");
  }

  const auto new_size = static_cast<size_t>(new_n) * static_cast<size_t>(new_n);
  std::vector<double> result(new_size, 0.0);

  for (int i = 0; i < old_n; ++i) {
    const auto src_start = m.begin() + static_cast<ptrdiff_t>(i * old_n);
    const auto src_end = src_start + old_n;
    auto dst_start = result.begin() + static_cast<ptrdiff_t>(i * new_n);
    std::copy(src_start, src_end, dst_start);
  }

  return result;
}

std::vector<double> MuhammadkhonIStressenAlgSEQ::UnpadMatrix(const std::vector<double> &m, int old_n, int new_n) {
  if (new_n > old_n) {
    throw std::runtime_error("new_n must be <= old_n for unpadding");
  }

  const auto new_size = static_cast<size_t>(new_n) * static_cast<size_t>(new_n);
  std::vector<double> result(new_size);

  for (int i = 0; i < new_n; ++i) {
    const auto src_start = m.begin() + static_cast<ptrdiff_t>(i * old_n);
    const auto src_end = src_start + new_n;
    auto dst_start = result.begin() + static_cast<ptrdiff_t>(i * new_n);
    std::copy(src_start, src_end, dst_start);
  }

  return result;
}

std::vector<double> MuhammadkhonIStressenAlgSEQ::Add(const std::vector<double> &a, const std::vector<double> &b,
                                                     int n) {
  const auto size = static_cast<size_t>(n) * static_cast<size_t>(n);
  if (a.size() != size || b.size() != size) {
    throw std::runtime_error("Invalid matrix sizes in Add");
  }

  std::vector<double> result(size);
  for (size_t i = 0; i < size; ++i) {
    result[i] = a[i] + b[i];
  }
  return result;
}

std::vector<double> MuhammadkhonIStressenAlgSEQ::Sub(const std::vector<double> &a, const std::vector<double> &b,
                                                     int n) {
  const auto size = static_cast<size_t>(n) * static_cast<size_t>(n);
  if (a.size() != size || b.size() != size) {
    throw std::runtime_error("Invalid matrix sizes in Sub");
  }

  std::vector<double> result(size);
  for (size_t i = 0; i < size; ++i) {
    result[i] = a[i] - b[i];
  }
  return result;
}

void MuhammadkhonIStressenAlgSEQ::Split(const std::vector<double> &parent, int n, std::vector<double> &a11,
                                        std::vector<double> &a12, std::vector<double> &a21, std::vector<double> &a22) {
  if (n % 2 != 0) {
    throw std::runtime_error("Matrix size must be even for Split");
  }

  const int half = n / 2;
  const auto half_size = static_cast<size_t>(half) * static_cast<size_t>(half);
  const auto parent_size = static_cast<size_t>(n) * static_cast<size_t>(n);

  if (parent.size() != parent_size) {
    throw std::runtime_error("Invalid parent matrix size in Split");
  }

  a11.resize(half_size);
  a12.resize(half_size);
  a21.resize(half_size);
  a22.resize(half_size);

  for (int i = 0; i < half; ++i) {
    for (int j = 0; j < half; ++j) {
      const auto idx = static_cast<size_t>((static_cast<ptrdiff_t>(i) * half) + j);
      const auto parent_idx_11 = static_cast<size_t>((static_cast<ptrdiff_t>(i) * n) + j);
      const auto parent_idx_12 = static_cast<size_t>((static_cast<ptrdiff_t>(i) * n) + j + half);
      const auto parent_idx_21 = static_cast<size_t>((static_cast<ptrdiff_t>(i + half) * n) + j);
      const auto parent_idx_22 = static_cast<size_t>((static_cast<ptrdiff_t>(i + half) * n) + j + half);

      a11[idx] = parent[parent_idx_11];
      a12[idx] = parent[parent_idx_12];
      a21[idx] = parent[parent_idx_21];
      a22[idx] = parent[parent_idx_22];
    }
  }
}

std::vector<double> MuhammadkhonIStressenAlgSEQ::Merge(const std::vector<double> &c11, const std::vector<double> &c12,
                                                       const std::vector<double> &c21, const std::vector<double> &c22,
                                                       int n) {
  const auto expected_size = static_cast<size_t>(n) * static_cast<size_t>(n);
  if (c11.size() != expected_size || c12.size() != expected_size || c21.size() != expected_size ||
      c22.size() != expected_size) {
    throw std::runtime_error("Invalid matrix sizes in Merge");
  }

  const int full = n * 2;
  const auto full_size = static_cast<size_t>(full) * static_cast<size_t>(full);
  std::vector<double> result(full_size);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      const auto src = static_cast<size_t>((static_cast<ptrdiff_t>(i) * n) + j);
      const auto dst_11 = static_cast<size_t>((static_cast<ptrdiff_t>(i) * full) + j);
      const auto dst_12 = static_cast<size_t>((static_cast<ptrdiff_t>(i) * full) + j + n);
      const auto dst_21 = static_cast<size_t>((static_cast<ptrdiff_t>(i + n) * full) + j);
      const auto dst_22 = static_cast<size_t>((static_cast<ptrdiff_t>(i + n) * full) + j + n);

      result[dst_11] = c11[src];
      result[dst_12] = c12[src];
      result[dst_21] = c21[src];
      result[dst_22] = c22[src];
    }
  }
  return result;
}

std::vector<double> MuhammadkhonIStressenAlgSEQ::NaiveMult(const std::vector<double> &a, const std::vector<double> &b,
                                                           int n) {
  const auto size = static_cast<size_t>(n) * static_cast<size_t>(n);
  const auto expected_size = static_cast<size_t>(n) * static_cast<size_t>(n);

  if (a.size() != expected_size || b.size() != expected_size) {
    throw std::runtime_error("Invalid matrix sizes in NaiveMult");
  }

  std::vector<double> c(size, 0.0);

  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < n; ++k) {
      const double aik = a[static_cast<size_t>((static_cast<ptrdiff_t>(i) * n) + k)];
      if (aik != 0.0) {
        for (int j = 0; j < n; ++j) {
          c[static_cast<size_t>((static_cast<ptrdiff_t>(i) * n) + j)] +=
              aik * b[static_cast<size_t>((static_cast<ptrdiff_t>(k) * n) + j)];
        }
      }
    }
  }
  return c;
}

std::vector<double> MuhammadkhonIStressenAlgSEQ::Strassen(const std::vector<double> &a, const std::vector<double> &b,
                                                          int n) {
  if (n <= kBaseCaseSize) {
    return NaiveMult(a, b, n);
  }

  if (n % 2 != 0) {
    throw std::runtime_error("Matrix size must be even for Strassen algorithm");
  }

  int half = n / 2;

  std::vector<double> a11, a12, a21, a22, b11, b12, b21, b22;
  Split(a, n, a11, a12, a21, a22);
  Split(b, n, b11, b12, b21, b22);

  auto p1 = Strassen(Add(a11, a22, half), Add(b11, b22, half), half);
  auto p2 = Strassen(Add(a21, a22, half), std::move(b11), half);
  auto p3 = Strassen(std::move(a11), Sub(b12, b22, half), half);
  auto p4 = Strassen(std::move(a22), Sub(b21, b11, half), half);
  auto p5 = Strassen(Add(a11, a12, half), std::move(b22), half);  // a11 уже перемещено в p3
  auto p6 = Strassen(Sub(a21, a11, half), Add(b11, b12, half), half);
  auto p7 = Strassen(Sub(a12, a22, half), Add(b21, b22, half), half);

  auto c11 = Add(Sub(Add(std::move(p1), std::move(p4), half), p5, half), p7, half);
  auto c12 = Add(std::move(p3), std::move(p5), half);
  auto c21 = Add(std::move(p2), std::move(p4), half);
  auto c22 = Add(Sub(Add(std::move(p1), std::move(p3), half), p2, half), p6, half);

  return Merge(std::move(c11), std::move(c12), std::move(c21), std::move(c22), half);
}

}  // namespace muhammadkhon_i_stressen_alg
