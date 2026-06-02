#include "akhmetov_daniil_strassen_dense_double/omp/include/ops_omp.hpp"

#include <cstddef>
#include <utility>
#include <vector>

#include "akhmetov_daniil_strassen_dense_double/common/include/common.hpp"

namespace akhmetov_daniil_strassen_dense_double {

AkhmetovDStrassenDenseDoubleOMP::AkhmetovDStrassenDenseDoubleOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool AkhmetovDStrassenDenseDoubleOMP::ValidationImpl() {
  const auto &input = GetInput();
  if (input.empty()) {
    return false;
  }
  const size_t n = format::GetN(input);
  if (n == 0) {
    return false;
  }
  const size_t expected_size = 1 + (2 * n * n);
  return input.size() == expected_size;
}

bool AkhmetovDStrassenDenseDoubleOMP::PreProcessingImpl() {
  const auto &input = GetInput();
  const size_t n = format::GetN(input);
  GetOutput().assign(n * n, 0.0);
  return true;
}

namespace {

constexpr size_t kThreshold = 64;
[[maybe_unused]] constexpr size_t kParallelThreshold = 256;

inline size_t NextPow2(size_t n) {
  size_t p = 1;
  while (p < n) {
    p <<= 1;
  }
  return p;
}

Matrix StandardMultiply(const Matrix &a, const Matrix &b, size_t size) {
  Matrix c(size * size, 0.0);
#pragma omp parallel for default(none) shared(a, b, c, size) schedule(static) if (size >= kParallelThreshold)
  for (size_t i = 0; i < size; ++i) {
    for (size_t k = 0; k < size; ++k) {
      const double aik = a.at((i * size) + k);
      const size_t bk = k * size;
      const size_t ci = i * size;
      for (size_t j = 0; j < size; ++j) {
        c.at(ci + j) += aik * b.at(bk + j);
      }
    }
  }
  return c;
}

void Split(const Matrix &src, Matrix &a11, Matrix &a12, Matrix &a21, Matrix &a22, size_t size) {
  const size_t half = size / 2;
  const size_t block = half * half;
  a11.assign(block, 0.0);
  a12.assign(block, 0.0);
  a21.assign(block, 0.0);
  a22.assign(block, 0.0);

#pragma omp parallel for default(none) shared(src, a11, a12, a21, a22, size, half) \
    schedule(static) if (size >= kParallelThreshold)
  for (size_t i = 0; i < half; ++i) {
    const size_t is = i * size;
    const size_t ih = i * half;
    const size_t is2 = (i + half) * size;
    for (size_t j = 0; j < half; ++j) {
      a11.at(ih + j) = src.at(is + j);
      a12.at(ih + j) = src.at(is + j + half);
      a21.at(ih + j) = src.at(is2 + j);
      a22.at(ih + j) = src.at(is2 + j + half);
    }
  }
}

void Merge(Matrix &dst, const Matrix &c11, const Matrix &c12, const Matrix &c21, const Matrix &c22, size_t size) {
  const size_t half = size / 2;
#pragma omp parallel for default(none) shared(dst, c11, c12, c21, c22, size, half) \
    schedule(static) if (size >= kParallelThreshold)
  for (size_t i = 0; i < half; ++i) {
    const size_t is = i * size;
    const size_t ih = i * half;
    const size_t is2 = (i + half) * size;
    for (size_t j = 0; j < half; ++j) {
      dst.at(is + j) = c11.at(ih + j);
      dst.at(is + j + half) = c12.at(ih + j);
      dst.at(is2 + j) = c21.at(ih + j);
      dst.at(is2 + j + half) = c22.at(ih + j);
    }
  }
}

inline void AddInto(const Matrix &a, const Matrix &b, Matrix &c) {
  const size_t n = a.size();
#pragma omp parallel for default(none) shared(a, b, c, n) \
    schedule(static) if (n >= (kParallelThreshold * kParallelThreshold))
  for (size_t i = 0; i < n; ++i) {
    c.at(i) = a.at(i) + b.at(i);
  }
}

inline void SubInto(const Matrix &a, const Matrix &b, Matrix &c) {
  const size_t n = a.size();
#pragma omp parallel for default(none) shared(a, b, c, n) \
    schedule(static) if (n >= (kParallelThreshold * kParallelThreshold))
  for (size_t i = 0; i < n; ++i) {
    c.at(i) = a.at(i) - b.at(i);
  }
}

struct Frame {
  Matrix a;
  Matrix b;
  Matrix result;
  size_t size;
  size_t stage{0};

  Matrix a11;
  Matrix a12;
  Matrix a21;
  Matrix a22;
  Matrix b11;
  Matrix b12;
  Matrix b21;
  Matrix b22;

  Matrix m1;
  Matrix m2;
  Matrix m3;
  Matrix m4;
  Matrix m5;
  Matrix m6;
  Matrix m7;

  Matrix temp_a;
  Matrix temp_b;

  Frame(Matrix aa, Matrix bb, size_t s) : a(std::move(aa)), b(std::move(bb)), size(s) {}
};

void ProcessTopFrame(std::vector<Frame> &stack, Matrix &final_result) {
  if (stack.empty()) {
    return;
  }

  const size_t current_index = stack.size() - 1;
  Frame &frame = stack.at(current_index);

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
      Split(frame.a, frame.a11, frame.a12, frame.a21, frame.a22, frame.size);
      Split(frame.b, frame.b11, frame.b12, frame.b21, frame.b22, frame.size);

      frame.temp_a.assign(block_size, 0.0);
      frame.temp_b.assign(block_size, 0.0);
      frame.m1.assign(block_size, 0.0);
      frame.m2.assign(block_size, 0.0);
      frame.m3.assign(block_size, 0.0);
      frame.m4.assign(block_size, 0.0);
      frame.m5.assign(block_size, 0.0);
      frame.m6.assign(block_size, 0.0);
      frame.m7.assign(block_size, 0.0);

      frame.stage = 1;
      return;
    }

    case 1: {
      AddInto(frame.a11, frame.a22, frame.temp_a);
      AddInto(frame.b11, frame.b22, frame.temp_b);

      stack.emplace_back(frame.temp_a, frame.temp_b, half);
      stack.at(current_index).stage = 2;
      return;
    }

    case 2: {
      frame.m1 = frame.temp_a;
      AddInto(frame.a21, frame.a22, frame.temp_a);

      stack.emplace_back(frame.temp_a, frame.b11, half);
      stack.at(current_index).stage = 3;
      return;
    }

    case 3: {
      frame.m2 = frame.temp_a;
      SubInto(frame.b12, frame.b22, frame.temp_b);

      stack.emplace_back(frame.a11, frame.temp_b, half);
      stack.at(current_index).stage = 4;
      return;
    }

    case 4: {
      frame.m3 = frame.temp_a;
      SubInto(frame.b21, frame.b11, frame.temp_b);

      stack.emplace_back(frame.a22, frame.temp_b, half);
      stack.at(current_index).stage = 5;
      return;
    }

    case 5: {
      frame.m4 = frame.temp_a;
      AddInto(frame.a11, frame.a12, frame.temp_a);

      stack.emplace_back(frame.temp_a, frame.b22, half);
      stack.at(current_index).stage = 6;
      return;
    }

    case 6: {
      frame.m5 = frame.temp_a;
      SubInto(frame.a21, frame.a11, frame.temp_a);
      AddInto(frame.b11, frame.b12, frame.temp_b);

      stack.emplace_back(frame.temp_a, frame.temp_b, half);
      stack.at(current_index).stage = 7;
      return;
    }

    case 7: {
      frame.m6 = frame.temp_a;
      SubInto(frame.a12, frame.a22, frame.temp_a);
      AddInto(frame.b21, frame.b22, frame.temp_b);

      stack.emplace_back(frame.temp_a, frame.temp_b, half);
      stack.at(current_index).stage = 8;
      return;
    }

    case 8: {
      frame.m7 = frame.temp_a;

      Matrix c11(block_size);
      Matrix c12(block_size);
      Matrix c21(block_size);
      Matrix c22(block_size);

#pragma omp parallel for default(none) shared(c11, c12, c21, c22, frame, block_size) \
    schedule(static) if (block_size >= (kParallelThreshold * kParallelThreshold))
      for (size_t i = 0; i < block_size; ++i) {
        c11.at(i) = frame.m1.at(i) + frame.m4.at(i) - frame.m5.at(i) + frame.m7.at(i);
        c12.at(i) = frame.m3.at(i) + frame.m5.at(i);
        c21.at(i) = frame.m2.at(i) + frame.m4.at(i);
        c22.at(i) = frame.m1.at(i) - frame.m2.at(i) + frame.m3.at(i) + frame.m6.at(i);
      }

      Matrix merged(frame.size * frame.size, 0.0);
      Merge(merged, c11, c12, c21, c22, frame.size);

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

Matrix PadTo(const Matrix &src, size_t n, size_t new_n) {
  if (new_n == n) {
    return src;
  }
  Matrix padded(new_n * new_n, 0.0);
#pragma omp parallel for default(none) shared(padded, src, n, new_n) schedule(static) if (new_n >= kParallelThreshold)
  for (size_t i = 0; i < n; ++i) {
    const size_t is = i * n;
    const size_t id = i * new_n;
    for (size_t j = 0; j < n; ++j) {
      padded.at(id + j) = src.at(is + j);
    }
  }
  return padded;
}

}  // namespace

bool AkhmetovDStrassenDenseDoubleOMP::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  const size_t n = format::GetN(input);
  const Matrix a = format::GetA(input);
  const Matrix b = format::GetB(input);

  if (n <= kThreshold) {
    output = StandardMultiply(a, b, n);
    return true;
  }

  const size_t new_n = NextPow2(n);
  const Matrix a_padded = PadTo(a, n, new_n);
  const Matrix b_padded = PadTo(b, n, new_n);

  const Matrix result_padded = StrassenMultiply(a_padded, b_padded, new_n);

  output.assign(n * n, 0.0);
#pragma omp parallel for default(none) shared(output, result_padded, n, new_n) \
    schedule(static) if (n >= kParallelThreshold)
  for (size_t i = 0; i < n; ++i) {
    const size_t is = i * new_n;
    const size_t id = i * n;
    for (size_t j = 0; j < n; ++j) {
      output.at(id + j) = result_padded.at(is + j);
    }
  }

  return true;
}

bool AkhmetovDStrassenDenseDoubleOMP::PostProcessingImpl() {
  const auto &input = GetInput();
  const size_t n = format::GetN(input);
  return GetOutput().size() == n * n;
}

}  // namespace akhmetov_daniil_strassen_dense_double
