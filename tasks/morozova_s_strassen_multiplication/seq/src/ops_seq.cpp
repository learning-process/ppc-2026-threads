#include "morozova_s_strassen_multiplication/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

#include "morozova_s_strassen_multiplication/common/include/common.hpp"

namespace morozova_s_strassen_multiplication {

namespace {

Matrix AddMatrixImpl(const Matrix &a, const Matrix &b) {
  int n = a.size;
  Matrix result(n);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      result(i, j) = a(i, j) + b(i, j);
    }
  }

  return result;
}

Matrix SubtractMatrixImpl(const Matrix &a, const Matrix &b) {
  int n = a.size;
  Matrix result(n);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      result(i, j) = a(i, j) - b(i, j);
    }
  }

  return result;
}

Matrix MultiplyStandardImpl(const Matrix &a, const Matrix &b) {
  int n = a.size;
  Matrix result(n);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double sum = 0.0;
      for (int k = 0; k < n; ++k) {
        sum += a(i, k) * b(k, j);
      }
      result(i, j) = sum;
    }
  }

  return result;
}

void SplitMatrixImpl(const Matrix &m, Matrix &m11, Matrix &m12, Matrix &m21, Matrix &m22) {
  int n = m.size;
  int half = n / 2;

  for (int i = 0; i < half; ++i) {
    for (int j = 0; j < half; ++j) {
      m11(i, j) = m(i, j);
      m12(i, j) = m(i, j + half);
      m21(i, j) = m(i + half, j);
      m22(i, j) = m(i + half, j + half);
    }
  }
}

Matrix MergeMatricesImpl(const Matrix &m11, const Matrix &m12, const Matrix &m21, const Matrix &m22) {
  int half = m11.size;
  int n = 2 * half;
  Matrix result(n);

  for (int i = 0; i < half; ++i) {
    for (int j = 0; j < half; ++j) {
      result(i, j) = m11(i, j);
      result(i, j + half) = m12(i, j);
      result(i + half, j) = m21(i, j);
      result(i + half, j + half) = m22(i, j);
    }
  }

  return result;
}

Matrix MultiplyStrassenIterative(const Matrix &a, const Matrix &b, int leaf_size) {
  int n = a.size;

  if (n <= leaf_size || n % 2 != 0) {
    return MultiplyStandardImpl(a, b);
  }

  int half = n / 2;

  Matrix a11(half);
  Matrix a12(half);
  Matrix a21(half);
  Matrix a22(half);
  Matrix b11(half);
  Matrix b12(half);
  Matrix b21(half);
  Matrix b22(half);

  SplitMatrixImpl(a, a11, a12, a21, a22);
  SplitMatrixImpl(b, b11, b12, b21, b22);

  Matrix p1 = MultiplyStandardImpl(a11, SubtractMatrixImpl(b12, b22));
  Matrix p2 = MultiplyStandardImpl(AddMatrixImpl(a11, a12), b22);
  Matrix p3 = MultiplyStandardImpl(AddMatrixImpl(a21, a22), b11);
  Matrix p4 = MultiplyStandardImpl(a22, SubtractMatrixImpl(b21, b11));
  Matrix p5 = MultiplyStandardImpl(AddMatrixImpl(a11, a22), AddMatrixImpl(b11, b22));
  Matrix p6 = MultiplyStandardImpl(SubtractMatrixImpl(a12, a22), AddMatrixImpl(b21, b22));
  Matrix p7 = MultiplyStandardImpl(SubtractMatrixImpl(a11, a21), AddMatrixImpl(b11, b12));

  Matrix c11 = AddMatrixImpl(SubtractMatrixImpl(AddMatrixImpl(p5, p4), p2), p6);
  Matrix c12 = AddMatrixImpl(p1, p2);
  Matrix c21 = AddMatrixImpl(p3, p4);
  Matrix c22 = SubtractMatrixImpl(SubtractMatrixImpl(AddMatrixImpl(p5, p1), p3), p7);

  return MergeMatricesImpl(c11, c12, c21, c22);
}

}  // namespace

MorozovaSStrassenMultiplicationSEQ::MorozovaSStrassenMultiplicationSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool MorozovaSStrassenMultiplicationSEQ::ValidationImpl() {
  return true;
}

bool MorozovaSStrassenMultiplicationSEQ::PreProcessingImpl() {
  if (GetInput().empty()) {
    valid_data_ = false;
    return true;
  }

  double size_val = GetInput()[0];
  if (size_val <= 0.0) {
    valid_data_ = false;
    return true;
  }

  int n = static_cast<int>(size_val);

  if (GetInput().size() != 1 + (2 * static_cast<size_t>(n) * static_cast<size_t>(n))) {
    valid_data_ = false;
    return true;
  }

  valid_data_ = true;
  n_ = n;

  a_ = Matrix(n_);
  b_ = Matrix(n_);

  int idx = 1;
  for (int i = 0; i < n_; ++i) {
    for (int j = 0; j < n_; ++j) {
      a_(i, j) = GetInput()[idx++];
    }
  }

  for (int i = 0; i < n_; ++i) {
    for (int j = 0; j < n_; ++j) {
      b_(i, j) = GetInput()[idx++];
    }
  }

  return true;
}

bool MorozovaSStrassenMultiplicationSEQ::RunImpl() {
  if (!valid_data_) {
    return true;
  }

  const int leaf_size = 64;

  if (n_ <= leaf_size) {
    c_ = MultiplyStandard(a_, b_);
  } else {
    c_ = MultiplyStrassenIterative(a_, b_, leaf_size);
  }

  return true;
}

bool MorozovaSStrassenMultiplicationSEQ::PostProcessingImpl() {
  OutType &output = GetOutput();
  output.clear();

  if (!valid_data_) {
    return true;
  }

  output.push_back(static_cast<double>(n_));

  for (int i = 0; i < n_; ++i) {
    for (int j = 0; j < n_; ++j) {
      output.push_back(c_(i, j));
    }
  }

  return true;
}

Matrix MorozovaSStrassenMultiplicationSEQ::AddMatrix(const Matrix &a, const Matrix &b) {
  return AddMatrixImpl(a, b);
}

Matrix MorozovaSStrassenMultiplicationSEQ::SubtractMatrix(const Matrix &a, const Matrix &b) {
  return SubtractMatrixImpl(a, b);
}

Matrix MorozovaSStrassenMultiplicationSEQ::MultiplyStandard(const Matrix &a, const Matrix &b) {
  return MultiplyStandardImpl(a, b);
}

void MorozovaSStrassenMultiplicationSEQ::SplitMatrix(const Matrix &m, Matrix &m11, Matrix &m12, Matrix &m21,
                                                     Matrix &m22) {
  SplitMatrixImpl(m, m11, m12, m21, m22);
}

Matrix MorozovaSStrassenMultiplicationSEQ::MergeMatrices(const Matrix &m11, const Matrix &m12, const Matrix &m21,
                                                         const Matrix &m22) {
  return MergeMatricesImpl(m11, m12, m21, m22);
}

Matrix MorozovaSStrassenMultiplicationSEQ::MultiplyStrassen(const Matrix &a, const Matrix &b, int leaf_size) {
  return MultiplyStrassenIterative(a, b, leaf_size);
}

}  // namespace morozova_s_strassen_multiplication
