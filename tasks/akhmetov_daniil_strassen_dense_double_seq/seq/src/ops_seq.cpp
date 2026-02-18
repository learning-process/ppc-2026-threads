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

  size_t half = size / 2;
  size_t half2 = half * half;

  Matrix a11(half2);
  Matrix a12(half2);
  Matrix a21(half2);
  Matrix a22(half2);

  Matrix b11(half2);
  Matrix b12(half2);
  Matrix b21(half2);
  Matrix b22(half2);

  SplitMatrix(a, a11, a12, a21, a22, size, half);
  SplitMatrix(b, b11, b12, b21, b22, size, half);

  Matrix temp1(half2);
  Matrix temp2(half2);

  for (size_t i = 0; i < half2; ++i) {
    temp1[i] = b12[i] - b22[i];
  }
  Matrix p1 = StrassenMultiply(a11, temp1, half);

  for (size_t i = 0; i < half2; ++i) {
    temp1[i] = a11[i] + a12[i];
  }
  Matrix p2 = StrassenMultiply(temp1, b22, half);

  for (size_t i = 0; i < half2; ++i) {
    temp1[i] = a21[i] + a22[i];
  }
  Matrix p3 = StrassenMultiply(temp1, b11, half);

  for (size_t i = 0; i < half2; ++i) {
    temp1[i] = b21[i] - b11[i];
  }
  Matrix p4 = StrassenMultiply(a22, temp1, half);

  for (size_t i = 0; i < half2; ++i) {
    temp1[i] = a11[i] + a22[i];
  }
  for (size_t i = 0; i < half2; ++i) {
    temp2[i] = b11[i] + b22[i];
  }
  Matrix p5 = StrassenMultiply(temp1, temp2, half);

  for (size_t i = 0; i < half2; ++i) {
    temp1[i] = a12[i] - a22[i];
  }
  for (size_t i = 0; i < half2; ++i) {
    temp2[i] = b21[i] + b22[i];
  }
  Matrix p6 = StrassenMultiply(temp1, temp2, half);

  for (size_t i = 0; i < half2; ++i) {
    temp1[i] = a11[i] - a21[i];
  }
  for (size_t i = 0; i < half2; ++i) {
    temp2[i] = b11[i] + b12[i];
  }
  Matrix p7 = StrassenMultiply(temp1, temp2, half);

  Matrix c11(half2);
  Matrix c12(half2);
  Matrix c21(half2);
  Matrix c22(half2);

  for (size_t i = 0; i < half2; ++i) {
    c11[i] = p5[i] + p4[i] - p2[i] + p6[i];
    c12[i] = p1[i] + p2[i];
    c21[i] = p3[i] + p4[i];
    c22[i] = p5[i] + p1[i] - p3[i] - p7[i];
  }

  Matrix c(size * size);
  MergeMatrix(c, c11, c12, c21, c22, size, half);
  return c;
}

}  // namespace

bool AkhmetovDStrassenDenseDoubleSEQ::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  size_t n = format::GetN(input);
  Matrix a = format::GetA(input);
  Matrix b = format::GetB(input);

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

  output.resize(n * n);
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
