#include "akhmetov_daniil_strassen_dense_double_seq/seq/include/ops_seq.hpp"

#include <cstddef>
#include <utility>
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

struct Frame {
  Matrix a;
  Matrix b;
  Matrix result;
  size_t size;
  size_t stage{0};

  Matrix a11, a12, a21, a22;
  Matrix b11, b12, b21, b22;

  Matrix m1, m2, m3, m4, m5, m6, m7;
  Matrix temp_a, temp_b;

  Frame(Matrix a, Matrix b, size_t s) : a(std::move(a)), b(std::move(b)), size(s) {}
};

void ProcessTopFrame(std::vector<Frame> &stack, Matrix &final_result) {
  if (stack.empty()) {
    return;
  }

  size_t current_index = stack.size() - 1;
  Frame &frame = stack[current_index];

  if (frame.size <= kThreshold) {
    Matrix base = StandardMultiply(frame.a, frame.b, frame.size);
    stack.pop_back();
    if (stack.empty()) {
      final_result = std::move(base);
    } else {
      stack.back().temp_a = std::move(base);
    }
    return;
  }

  const size_t half = frame.size / 2;
  const size_t block_size = half * half;

  switch (frame.stage) {
    case 0: {
      frame.a11.resize(block_size);
      frame.a12.resize(block_size);
      frame.a21.resize(block_size);
      frame.a22.resize(block_size);
      frame.b11.resize(block_size);
      frame.b12.resize(block_size);
      frame.b21.resize(block_size);
      frame.b22.resize(block_size);

      SplitMatrix(frame.a, frame.a11, frame.a12, frame.a21, frame.a22, frame.size, half);
      SplitMatrix(frame.b, frame.b11, frame.b12, frame.b21, frame.b22, frame.size, half);

      frame.temp_a.resize(block_size);
      frame.temp_b.resize(block_size);
      frame.m1.resize(block_size);
      frame.m2.resize(block_size);
      frame.m3.resize(block_size);
      frame.m4.resize(block_size);
      frame.m5.resize(block_size);
      frame.m6.resize(block_size);
      frame.m7.resize(block_size);

      frame.stage = 1;
      return;
    }

    case 1: {
      Add(frame.a11, frame.a22, frame.temp_a);
      Add(frame.b11, frame.b22, frame.temp_b);

      Matrix temp_a_copy = frame.temp_a;
      Matrix temp_b_copy = frame.temp_b;

      stack.emplace_back(std::move(temp_a_copy), std::move(temp_b_copy), half);

      Frame &updated_frame = stack[current_index];
      updated_frame.stage = 2;
      return;
    }

    case 2: {
      frame.m1 = frame.temp_a;
      Add(frame.a21, frame.a22, frame.temp_a);

      Matrix temp_a_copy = frame.temp_a;
      stack.emplace_back(std::move(temp_a_copy), frame.b11, half);

      Frame &updated_frame = stack[current_index];
      updated_frame.stage = 3;
      return;
    }

    case 3: {
      frame.m2 = frame.temp_a;
      Sub(frame.b12, frame.b22, frame.temp_b);

      Matrix temp_b_copy = frame.temp_b;
      stack.emplace_back(frame.a11, std::move(temp_b_copy), half);

      Frame &updated_frame = stack[current_index];
      updated_frame.stage = 4;
      return;
    }

    case 4: {
      frame.m3 = frame.temp_a;
      Sub(frame.b21, frame.b11, frame.temp_b);

      Matrix temp_b_copy = frame.temp_b;
      stack.emplace_back(frame.a22, std::move(temp_b_copy), half);

      Frame &updated_frame = stack[current_index];
      updated_frame.stage = 5;
      return;
    }

    case 5: {
      frame.m4 = frame.temp_a;
      Add(frame.a11, frame.a12, frame.temp_a);

      Matrix temp_a_copy = frame.temp_a;
      stack.emplace_back(std::move(temp_a_copy), frame.b22, half);

      Frame &updated_frame = stack[current_index];
      updated_frame.stage = 6;
      return;
    }

    case 6: {
      frame.m5 = frame.temp_a;
      Sub(frame.a21, frame.a11, frame.temp_a);
      Add(frame.b11, frame.b12, frame.temp_b);

      Matrix temp_a_copy = frame.temp_a;
      Matrix temp_b_copy = frame.temp_b;
      stack.emplace_back(std::move(temp_a_copy), std::move(temp_b_copy), half);

      Frame &updated_frame = stack[current_index];
      updated_frame.stage = 7;
      return;
    }

    case 7: {
      frame.m6 = frame.temp_a;
      Sub(frame.a12, frame.a22, frame.temp_a);
      Add(frame.b21, frame.b22, frame.temp_b);

      Matrix temp_a_copy = frame.temp_a;
      Matrix temp_b_copy = frame.temp_b;
      stack.emplace_back(std::move(temp_a_copy), std::move(temp_b_copy), half);

      Frame &updated_frame = stack[current_index];
      updated_frame.stage = 8;
      return;
    }

    case 8: {
      frame.m7 = frame.temp_a;

      Matrix c11(block_size);
      Matrix c12(block_size);
      Matrix c21(block_size);
      Matrix c22(block_size);

      for (size_t i = 0; i < block_size; ++i) {
        c11[i] = frame.m1[i] + frame.m4[i] - frame.m5[i] + frame.m7[i];
        c12[i] = frame.m3[i] + frame.m5[i];
        c21[i] = frame.m2[i] + frame.m4[i];
        c22[i] = frame.m1[i] - frame.m2[i] + frame.m3[i] + frame.m6[i];
      }

      Matrix merged(frame.size * frame.size);
      MergeMatrix(merged, c11, c12, c21, c22, frame.size, half);

      stack.pop_back();
      if (stack.empty()) {
        final_result = std::move(merged);
      } else {
        stack.back().temp_a = std::move(merged);
      }
      return;
    }

    default:
      return;
  }
}

Matrix StrassenMultiply(const Matrix &a_init, const Matrix &b_init, size_t size_init) {
  std::vector<Frame> stack;
  stack.emplace_back(a_init, b_init, size_init);

  Matrix result;

  while (!stack.empty()) {
    ProcessTopFrame(stack, result);
  }

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
