#include "akhmetov_daniil_strassen_dense_double_seq/seq/include/ops_seq.hpp"

#include <cstddef>
#include <vector>

#include "akhmetov_daniil_strassen_dense_double_seq/common/include/common.hpp"

namespace akhmetov_daniil_strassen_dense_double_seq {

AkhmetovDStrassenDenseDoubleSEQ::AkhmetovDStrassenDenseDoubleSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool AkhmetovDStrassenDenseDoubleSEQ::ValidationImpl() {
  const auto &input = GetInput();
  if (input.empty()) {
    return false;
  }
  size_t n = format::GetN(input);
  if (n == 0) {
    return false;
  }
  size_t expected_size = 1 + (2 * n * n);
  return input.size() == expected_size;
}

bool AkhmetovDStrassenDenseDoubleSEQ::PreProcessingImpl() {
  const auto &input = GetInput();
  size_t n = format::GetN(input);
  GetOutput().resize(n * n);
  return true;
}

namespace {
const size_t kThreshold = 64;

Matrix StandardMultiply(const Matrix &a, const Matrix &b, size_t size) {
  Matrix c(size * size, 0.0);
  for (size_t i = 0; i < size; ++i) {
    for (size_t j = 0; j < size; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < size; ++k) {
        sum += a[(i * size) + k] * b[(k * size) + j];
      }
      c[(i * size) + j] = sum;
    }
  }
  return c;
}

inline void Add(const Matrix &a, const Matrix &b, Matrix &c) {
  const size_t n = a.size();
  for (size_t i = 0; i < n; ++i) {
    c[i] = a[i] + b[i];
  }
}

inline void Sub(const Matrix &a, const Matrix &b, Matrix &c) {
  const size_t n = a.size();
  for (size_t i = 0; i < n; ++i) {
    c[i] = a[i] - b[i];
  }
}

void SplitMatrix(const Matrix &src, Matrix &a11, Matrix &a12, Matrix &a21, Matrix &a22, size_t size, size_t half) {
  for (size_t i = 0; i < half; ++i) {
    for (size_t j = 0; j < half; ++j) {
      a11[(i * half) + j] = src[(i * size) + j];
      a12[(i * half) + j] = src[(i * size) + j + half];
      a21[(i * half) + j] = src[((i + half) * size) + j];
      a22[(i * half) + j] = src[((i + half) * size) + j + half];
    }
  }
}

void MergeMatrix(Matrix &dst, const Matrix &c11, const Matrix &c12, const Matrix &c21, const Matrix &c22, size_t size,
                 size_t half) {
  for (size_t i = 0; i < half; ++i) {
    for (size_t j = 0; j < half; ++j) {
      dst[(i * size) + j] = c11[(i * half) + j];
      dst[(i * size) + j + half] = c12[(i * half) + j];
      dst[((i + half) * size) + j] = c21[(i * half) + j];
      dst[((i + half) * size) + j + half] = c22[(i * half) + j];
    }
  }
}

void CreatePaddedMatrices(const Matrix &a, const Matrix &b, size_t n, size_t new_n, Matrix &a_padded,
                          Matrix &b_padded) {
  if (new_n != n) {
    a_padded.assign(new_n * new_n, 0.0);
    b_padded.assign(new_n * new_n, 0.0);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        a_padded[(i * new_n) + j] = a[(i * n) + j];
        b_padded[(i * new_n) + j] = b[(i * n) + j];
      }
    }
  } else {
    a_padded = a;
    b_padded = b;
  }
}

// NOLINTNEXTLINE(misc-no-recursion)
Matrix StrassenMultiply(const Matrix &a, const Matrix &b, size_t size) {
  if (size <= kThreshold) {
    return StandardMultiply(a, b, size);
  }

  const size_t half = size / 2;
  const size_t block_size = half * half;

  Matrix a11(block_size);
  Matrix a12(block_size);
  Matrix a21(block_size);
  Matrix a22(block_size);
  Matrix b11(block_size);
  Matrix b12(block_size);
  Matrix b21(block_size);
  Matrix b22(block_size);

  SplitMatrix(a, a11, a12, a21, a22, size, half);
  SplitMatrix(b, b11, b12, b21, b22, size, half);

  Matrix temp_a(block_size);
  Matrix temp_b(block_size);

  // M1 = (A11 + A22)(B11 + B22)
  Add(a11, a22, temp_a);
  Add(b11, b22, temp_b);
  Matrix m1 = StrassenMultiply(temp_a, temp_b, half);

  // M2 = (A21 + A22)B11
  Add(a21, a22, temp_a);
  Matrix m2 = StrassenMultiply(temp_a, b11, half);

  // M3 = A11(B12 - B22)
  Sub(b12, b22, temp_b);
  Matrix m3 = StrassenMultiply(a11, temp_b, half);

  // M4 = A22(B21 - B11)
  Sub(b21, b11, temp_b);
  Matrix m4 = StrassenMultiply(a22, temp_b, half);

  // M5 = (A11 + A12)B22
  Add(a11, a12, temp_a);
  Matrix m5 = StrassenMultiply(temp_a, b22, half);

  // M6 = (A21 - A11)(B11 + B12)
  Sub(a21, a11, temp_a);
  Add(b11, b12, temp_b);
  Matrix m6 = StrassenMultiply(temp_a, temp_b, half);

  // M7 = (A12 - A22)(B21 + B22)
  Sub(a12, a22, temp_a);
  Add(b21, b22, temp_b);
  Matrix m7 = StrassenMultiply(temp_a, temp_b, half);

  Matrix c11(block_size);
  Matrix c12(block_size);
  Matrix c21(block_size);
  Matrix c22(block_size);

  for (size_t i = 0; i < block_size; ++i) {
    c11[i] = m1[i] + m4[i] - m5[i] + m7[i];
    c12[i] = m3[i] + m5[i];
    c21[i] = m2[i] + m4[i];
    c22[i] = m1[i] - m2[i] + m3[i] + m6[i];
  }

  Matrix result(size * size);
  MergeMatrix(result, c11, c12, c21, c22, size, half);
  return result;
}

}  // namespace

bool AkhmetovDStrassenDenseDoubleSEQ::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  const size_t n = format::GetN(input);
  const Matrix a = format::GetA(input);
  const Matrix b = format::GetB(input);

  if (n <= kThreshold) {
    output = StandardMultiply(a, b, n);
    return true;
  }

  size_t new_n = 1;
  while (new_n < n) {
    new_n <<= 1;
  }

  Matrix a_padded;
  Matrix b_padded;
  CreatePaddedMatrices(a, b, n, new_n, a_padded, b_padded);

  Matrix result_padded = StrassenMultiply(a_padded, b_padded, new_n);

  output.assign(n * n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      output[(i * n) + j] = result_padded[(i * new_n) + j];
    }
  }

  return true;
}

bool AkhmetovDStrassenDenseDoubleSEQ::PostProcessingImpl() {
  const auto &input = GetInput();
  size_t n = format::GetN(input);
  return GetOutput().size() == n * n;
}

}  // namespace akhmetov_daniil_strassen_dense_double_seq
