#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/common/include/common.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/omp/include/ops_omp.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/seq/include/ops_seq.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {
namespace {
constexpr double kTestTol = 1e-12;

bool ComplexNear(const Complex &lhs, const Complex &rhs, double tol = kTestTol) {
  return std::abs(lhs.real() - rhs.real()) <= tol && std::abs(lhs.imag() - rhs.imag()) <= tol;
}

bool CrsEqual(const CRSMatrix &expected, const CRSMatrix &actual, double tol = kTestTol) {
  if (expected.rows != actual.rows || expected.cols != actual.cols) {
    return false;
  }
  if (expected.row_ptr != actual.row_ptr || expected.col_indices != actual.col_indices) {
    return false;
  }
  if (expected.values.size() != actual.values.size()) {
    return false;
  }
  for (std::size_t i = 0; i < expected.values.size(); ++i) {
    if (!ComplexNear(expected.values[i], actual.values[i], tol)) {
      return false;
    }
  }
  return true;
}

std::vector<std::vector<Complex>> ToDense(const CRSMatrix &m) {
  std::vector<std::vector<Complex>> dense(m.rows, std::vector<Complex>(m.cols, {0.0, 0.0}));
  for (int i = 0; i < m.rows; ++i) {
    for (int j = m.row_ptr[i]; j < m.row_ptr[i + 1]; ++j) {
      dense[i][m.col_indices[j]] = m.values[j];
    }
  }
  return dense;
}

bool DenseEqual(const std::vector<std::vector<Complex>> &expected, const CRSMatrix &actual, double tol = kTestTol) {
  const int rows = static_cast<int>(expected.size());
  const int cols = rows > 0 ? static_cast<int>(expected[0].size()) : 0;
  if (actual.rows != rows || actual.cols != cols) {
    return false;
  }
  const auto dense = ToDense(actual);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (!ComplexNear(dense[i][j], expected[i][j], tol)) {
        return false;
      }
    }
  }
  return true;
}

bool Dense2x2Equal(const CRSMatrix &lhs, const CRSMatrix &rhs, double tol = 1e-11) {
  const auto l = ToDense(lhs);
  const auto r = ToDense(rhs);
  if (l.size() != 2 || r.size() != 2 || l[0].size() != 2 || r[0].size() != 2) {
    return false;
  }
  return ComplexNear(l[0][0], r[0][0], tol) && ComplexNear(l[0][1], r[0][1], tol) &&
         ComplexNear(l[1][0], r[1][0], tol) && ComplexNear(l[1][1], r[1][1], tol);
}

bool CheckPartialCancellation(const CRSMatrix &c) {
  return c.rows == 2 && c.cols == 1 && c.row_ptr == std::vector<int>({0, 0, 1}) && c.values.size() == 1 &&
         ComplexNear(c.values[0], Complex(0.0, 2.0));
}

bool CheckCornerElements5x5(const CRSMatrix &c) {
  const auto dense = ToDense(c);
  return c.rows == 5 && c.cols == 5 && c.NonZeros() == 2U && ComplexNear(dense[0][0], Complex(2.0, 0.0)) &&
         ComplexNear(dense[4][4], Complex(1.0, 0.0));
}

bool CheckDenseRowTimesIdentity(const CRSMatrix &c) {
  const auto dense = ToDense(c);
  return c.NonZeros() == 4U && ComplexNear(dense[0][0], Complex(1.0, 1.0)) &&
         ComplexNear(dense[0][1], Complex(2.0, 0.0)) && ComplexNear(dense[0][2], Complex(3.0, -1.0)) &&
         ComplexNear(dense[0][3], Complex(0.0, 4.0));
}

bool IsShapeAndEmpty(const CRSMatrix &c, int rows, int cols) {
  return c.rows == rows && c.cols == cols && c.values.empty();
}

bool IsShapeAndValid(const CRSMatrix &c, int rows, int cols) {
  return c.rows == rows && c.cols == cols && c.IsValid();
}

bool IsSortedInRows(const CRSMatrix &c) {
  for (int i = 0; i < c.rows; ++i) {
    for (int j = c.row_ptr[i]; j < c.row_ptr[i + 1] - 1; ++j) {
      if (c.col_indices[j] >= c.col_indices[j + 1]) {
        return false;
      }
    }
  }
  return true;
}

bool HasNoExplicitZeros(const CRSMatrix &c) {
  return std::ranges::all_of(c.values, [](const Complex &v) { return std::abs(v) > kEpsilon; });
}

CRSMatrix RunSeqMultiply(const CRSMatrix &a, const CRSMatrix &b) {
  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(a, b));
  if (!(task.Validation() && task.PreProcessing() && task.Run() && task.PostProcessing())) {
    return CRSMatrix{};
  }
  return task.GetOutput();
}

struct TestCase {
  std::string name;
  CRSMatrix a;
  CRSMatrix b;
  std::function<bool(const CRSMatrix &)> check;
  bool expect_valid = true;
};

using TestCaseType = TestCase;

TestCaseType MakeInvalid(const std::string &name, CRSMatrix a, CRSMatrix b) {
  return {name, std::move(a), std::move(b), {}, false};
}

TestCaseType MakeCase(const std::string &name, CRSMatrix a, CRSMatrix b, std::function<bool(const CRSMatrix &)> check) {
  return {name, std::move(a), std::move(b), std::move(check), true};
}

const std::vector<std::vector<Complex>> kDense2x2Expected = {{Complex(3, 3), Complex(-1, 1)},
                                                             {Complex(6, 0), Complex(0, 2)}};
}  // namespace

class VidermanRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestCaseType> {
 public:
  static std::string PrintTestParam(const TestCaseType &test_param) {
    return test_param.name;
  }

 protected:
  void SetUp() override {
    const TestCaseType &params =
        std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    test_case_ = params;
    input_data_ = std::make_tuple(test_case_.a, test_case_.b);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (!test_case_.check) {
      return false;
    }
    return test_case_.check(output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  TestCaseType test_case_;
};

namespace {

TEST_P(VidermanRunFuncTests, CRSComplexMult) {
  const auto &test_param = GetParam();
  const auto &test_case = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(test_param);

  if (!test_case.expect_valid) {
    auto task =
        std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTaskGetter)>(test_param)(GetTestInputData());
    EXPECT_FALSE(task->Validation());
    return;
  }

  ExecuteTest(test_param);
}

const std::array<TestCaseType, 34> kTestParam = {
    // Validation (invalid cases)
    MakeInvalid("incompatible_dimensions", CRSMatrix(2, 3), CRSMatrix(4, 5)),
    [] {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, -1};
  a.values = {Complex(1, 0), Complex(2, 0)};

  CRSMatrix b(2, 2);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 1};
  b.values = {Complex(1, 0), Complex(1, 0)};
  return MakeInvalid("negative_col_index", a, b);
}(),
    [] {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 5};
  a.values = {Complex(1, 0), Complex(2, 0)};

  CRSMatrix b(2, 3);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 1};
  b.values = {Complex(1, 0), Complex(1, 0)};
  return MakeInvalid("col_index_out_of_range", a, b);
}(),
    [] {
  CRSMatrix a(1, 3);
  a.row_ptr = {0, 3};
  a.col_indices = {2, 0, 1};
  a.values = {Complex(1, 0), Complex(2, 0), Complex(3, 0)};

  CRSMatrix b(3, 3);
  b.row_ptr = {0, 1, 2, 3};
  b.col_indices = {0, 1, 2};
  b.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};
  return MakeInvalid("unsorted_col_indices", a, b);
}(),
    [] {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 2, 1};
  a.col_indices = {0, 1, 0};
  a.values = {Complex(1, 0), Complex(2, 0), Complex(3, 0)};

  CRSMatrix b(2, 2);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 1};
  b.values = {Complex(1, 0), Complex(1, 0)};
  return MakeInvalid("non_monotonic_row_ptr", a, b);
}(),
    [] {
  CRSMatrix a;
  a.rows = 3;
  a.cols = 3;
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 1};
  a.values = {Complex(1, 0), Complex(1, 0)};

  CRSMatrix b(3, 3);
  b.row_ptr = {0, 1, 2, 3};
  b.col_indices = {0, 1, 2};
  b.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};
  return MakeInvalid("wrong_row_ptr_size", a, b);
}(),
    [] {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 1};
  a.values = {Complex(1, 0)};

  CRSMatrix b(2, 2);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 1};
  b.values = {Complex(1, 0), Complex(1, 0)};
  return MakeInvalid("col_indices_values_size_mismatch", a, b);
}(),

    // Edge cases
    [] {
  CRSMatrix a(1, 1);
  a.row_ptr = {0, 1};
  a.col_indices = {0};
  a.values = {Complex(3.0, 4.0)};

  CRSMatrix b(1, 1);
  b.row_ptr = {0, 1};
  b.col_indices = {0};
  b.values = {Complex(1.0, -2.0)};

  CRSMatrix expected(1, 1);
  expected.row_ptr = {0, 1};
  expected.col_indices = {0};
  expected.values = {Complex(11.0, -2.0)};

  return MakeCase("single_element", a, b, [expected](const CRSMatrix &c) { return CrsEqual(expected, c); });
}(),
    MakeCase("both_zero_matrices", CRSMatrix(3, 4), CRSMatrix(4, 5),
             [](const CRSMatrix &c) { return IsShapeAndValid(c, 3, 5); }),
    [] {
  CRSMatrix a(2, 3);
  CRSMatrix b(3, 2);
  b.row_ptr = {0, 1, 2, 3};
  b.col_indices = {0, 1, 0};
  b.values = {Complex(1, 0), Complex(2, 0), Complex(3, 0)};
  return MakeCase("zero_a_nonzero_b", a, b, [](const CRSMatrix &c) { return IsShapeAndEmpty(c, 2, 2); });
}(),
    [] {
  CRSMatrix a(2, 3);
  a.row_ptr = {0, 2, 3};
  a.col_indices = {0, 2, 1};
  a.values = {Complex(1, 1), Complex(2, 0), Complex(3, -1)};
  return MakeCase("nonzero_a_zero_b", a, CRSMatrix(3, 4), [](const CRSMatrix &c) { return IsShapeAndEmpty(c, 2, 4); });
}(),
    [] {
  CRSMatrix a(1, 3);
  a.row_ptr = {0, 3};
  a.col_indices = {0, 1, 2};
  a.values = {Complex(1, 1), Complex(2, 0), Complex(3, -1)};

  CRSMatrix b(3, 1);
  b.row_ptr = {0, 1, 2, 3};
  b.col_indices = {0, 0, 0};
  b.values = {Complex(1, 0), Complex(0, 1), Complex(1, 1)};

  CRSMatrix expected(1, 1);
  expected.row_ptr = {0, 1};
  expected.col_indices = {0};
  expected.values = {Complex(5.0, 5.0)};
  return MakeCase("row_vector_times_col_vector", a, b,
                  [expected](const CRSMatrix &c) { return CrsEqual(expected, c); });
}(),
    [] {
  CRSMatrix a(2, 1);
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 0};
  a.values = {Complex(1, 1), Complex(2, 0)};

  CRSMatrix b(1, 2);
  b.row_ptr = {0, 2};
  b.col_indices = {0, 1};
  b.values = {Complex(3, 0), Complex(0, 1)};

  return MakeCase("col_vector_times_row_vector", a, b,
                  [](const CRSMatrix &c) { return DenseEqual(kDense2x2Expected, c); });
}(),
    [] {
  CRSMatrix a(4, 1);
  a.row_ptr = {0, 1, 2, 3, 4};
  a.col_indices = {0, 0, 0, 0};
  a.values = {Complex(1, 0), Complex(0, 1), Complex(-1, 0), Complex(0, -1)};

  CRSMatrix b(1, 4);
  b.row_ptr = {0, 4};
  b.col_indices = {0, 1, 2, 3};
  b.values = {Complex(1, 0), Complex(0, 1), Complex(-1, 0), Complex(0, -1)};

  std::vector<std::vector<Complex>> expected(4, std::vector<Complex>(4));
  const std::vector<Complex> a_vals = {Complex(1, 0), Complex(0, 1), Complex(-1, 0), Complex(0, -1)};
  const std::vector<Complex> b_vals = {Complex(1, 0), Complex(0, 1), Complex(-1, 0), Complex(0, -1)};
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      expected[i][j] = a_vals[i] * b_vals[j];
    }
  }
  return MakeCase("tall_skinny_matrix", a, b, [expected](const CRSMatrix &c) { return DenseEqual(expected, c); });
}(),

    // Complex arithmetic
    [] {
  CRSMatrix a(3, 3);
  a.row_ptr = {0, 1, 2, 3};
  a.col_indices = {0, 1, 2};
  a.values = {Complex(1, 1), Complex(2, 2), Complex(3, 3)};

  CRSMatrix b(3, 3);
  b.row_ptr = {0, 1, 2, 3};
  b.col_indices = {0, 1, 2};
  b.values = {Complex(4, 0), Complex(5, 0), Complex(6, 0)};

  CRSMatrix expected(3, 3);
  expected.row_ptr = {0, 1, 2, 3};
  expected.col_indices = {0, 1, 2};
  expected.values = {Complex(4, 4), Complex(10, 10), Complex(18, 18)};
  return MakeCase("diagonal_multiplication", a, b, [expected](const CRSMatrix &c) { return CrsEqual(expected, c); });
}(),
    [] {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 1};
  a.values = {Complex(0, 1), Complex(0, 1)};

  CRSMatrix b(2, 2);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 1};
  b.values = {Complex(0, 1), Complex(0, 1)};

  CRSMatrix expected(2, 2);
  expected.row_ptr = {0, 1, 2};
  expected.col_indices = {0, 1};
  expected.values = {Complex(-1, 0), Complex(-1, 0)};
  return MakeCase("pure_imaginary_squared", a, b, [expected](const CRSMatrix &c) { return CrsEqual(expected, c); });
}(),
    [] {
  CRSMatrix a(1, 1);
  a.row_ptr = {0, 1};
  a.col_indices = {0};
  a.values = {Complex(3, 4)};

  CRSMatrix b(1, 1);
  b.row_ptr = {0, 1};
  b.col_indices = {0};
  b.values = {Complex(3, -4)};

  CRSMatrix expected(1, 1);
  expected.row_ptr = {0, 1};
  expected.col_indices = {0};
  expected.values = {Complex(25, 0)};
  return MakeCase("conjugate_product", a, b, [expected](const CRSMatrix &c) { return CrsEqual(expected, c); });
}(),
    [] {
  CRSMatrix a(1, 2);
  a.row_ptr = {0, 2};
  a.col_indices = {0, 1};
  a.values = {Complex(1, 0), Complex(-1, 0)};

  CRSMatrix b(2, 1);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 0};
  b.values = {Complex(0, 1), Complex(0, 1)};
  return MakeCase("cancellation_to_zero", a, b, [](const CRSMatrix &c) { return IsShapeAndEmpty(c, 1, 1); });
}(),
    [] {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 2, 4};
  a.col_indices = {0, 1, 0, 1};
  a.values = {Complex(1, 0), Complex(-1, 0), Complex(1, 0), Complex(1, 0)};

  CRSMatrix b(2, 1);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 0};
  b.values = {Complex(0, 1), Complex(0, 1)};
  return MakeCase("partial_cancellation", a, b, [](const CRSMatrix &c) { return CheckPartialCancellation(c); });
}(),

    // Algebraic
    [] {
  CRSMatrix a(3, 3);
  a.row_ptr = {0, 2, 3, 4};
  a.col_indices = {0, 2, 1, 0};
  a.values = {Complex(1, 2), Complex(3, 0), Complex(0, 1), Complex(5, -1)};

  CRSMatrix i(3, 3);
  i.row_ptr = {0, 1, 2, 3};
  i.col_indices = {0, 1, 2};
  i.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};
  return MakeCase("right_identity", a, i, [a](const CRSMatrix &c) { return CrsEqual(a, c); });
}(),
    [] {
  CRSMatrix a(3, 3);
  a.row_ptr = {0, 2, 3, 4};
  a.col_indices = {0, 2, 1, 0};
  a.values = {Complex(1, 2), Complex(3, 0), Complex(0, 1), Complex(5, -1)};

  CRSMatrix i(3, 3);
  i.row_ptr = {0, 1, 2, 3};
  i.col_indices = {0, 1, 2};
  i.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};
  return MakeCase("left_identity", i, a, [a](const CRSMatrix &c) { return CrsEqual(a, c); });
}(),
    [] {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 2, 3};
  a.col_indices = {0, 1, 0};
  a.values = {Complex(1, 1), Complex(2, 0), Complex(0, 1)};

  CRSMatrix b(2, 2);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 1};
  b.values = {Complex(3, 0), Complex(0, 2)};

  CRSMatrix c(2, 2);
  c.row_ptr = {0, 2, 2};
  c.col_indices = {0, 1};
  c.values = {Complex(1, 0), Complex(1, 1)};

  return MakeCase("associativity", a, b, [a, b, c](const CRSMatrix &ab) {
    const CRSMatrix left = RunSeqMultiply(ab, c);
    const CRSMatrix right = RunSeqMultiply(a, RunSeqMultiply(b, c));
    return Dense2x2Equal(left, right);
  });
}(),
    [] {
  CRSMatrix i_i(2, 2);
  i_i.row_ptr = {0, 1, 2};
  i_i.col_indices = {0, 1};
  i_i.values = {Complex(0, 1), Complex(0, 1)};

  CRSMatrix minus_i(2, 2);
  minus_i.row_ptr = {0, 1, 2};
  minus_i.col_indices = {0, 1};
  minus_i.values = {Complex(-1, 0), Complex(-1, 0)};
  return MakeCase("square_of_scaled_identity", i_i, i_i,
                  [minus_i](const CRSMatrix &c) { return CrsEqual(minus_i, c); });
}(),
    [] {
  CRSMatrix p(3, 3);
  p.row_ptr = {0, 1, 2, 3};
  p.col_indices = {1, 2, 0};
  p.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  CRSMatrix pt(3, 3);
  pt.row_ptr = {0, 1, 2, 3};
  pt.col_indices = {2, 0, 1};
  pt.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  CRSMatrix i(3, 3);
  i.row_ptr = {0, 1, 2, 3};
  i.col_indices = {0, 1, 2};
  i.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};
  return MakeCase("permutation_times_transpose", p, pt, [i](const CRSMatrix &c) { return CrsEqual(i, c); });
}(),

    // Structural
    [] {
  CRSMatrix a(4, 3);
  a.row_ptr = {0, 2, 3, 3, 4};
  a.col_indices = {0, 2, 1, 2};
  a.values = {Complex(1, 1), Complex(2, 0), Complex(0, 3), Complex(-1, 1)};

  CRSMatrix b(3, 5);
  b.row_ptr = {0, 2, 3, 5};
  b.col_indices = {0, 3, 2, 1, 4};
  b.values = {Complex(1, 0), Complex(2, 1), Complex(0, 1), Complex(3, 0), Complex(1, -1)};
  return MakeCase("output_is_valid_crs", a, b, [](const CRSMatrix &c) { return IsShapeAndValid(c, 4, 5); });
}(),
    [] {
  CRSMatrix a(2, 3);
  a.row_ptr = {0, 3, 3};
  a.col_indices = {0, 1, 2};
  a.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  CRSMatrix b(3, 3);
  b.row_ptr = {0, 1, 2, 3};
  b.col_indices = {2, 0, 1};
  b.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};
  return MakeCase("col_indices_sorted_in_output", a, b, [](const CRSMatrix &c) { return IsSortedInRows(c); });
}(),
    [] {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 1};
  a.values = {Complex(1, 0), Complex(2, 0)};

  CRSMatrix b(2, 2);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 1};
  b.values = {Complex(3, 0), Complex(4, 0)};
  return MakeCase("row_ptr_starts_at_zero", a, b,
                  [](const CRSMatrix &c) { return !c.row_ptr.empty() && c.row_ptr[0] == 0; });
}(),
    [] {
  CRSMatrix a(3, 2);
  a.row_ptr = {0, 2, 2, 2};
  a.col_indices = {0, 1};
  a.values = {Complex(1, 1), Complex(2, 0)};

  CRSMatrix b(2, 3);
  b.row_ptr = {0, 2, 3};
  b.col_indices = {0, 2, 1};
  b.values = {Complex(1, 0), Complex(0, 1), Complex(3, 0)};
  return MakeCase("row_ptr_last_equals_nnz", a, b,
                  [](const CRSMatrix &c) { return static_cast<std::size_t>(c.row_ptr[c.rows]) == c.values.size(); });
}(),
    [] {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 1};
  a.values = {Complex(1, 0), Complex(1, 0)};

  CRSMatrix b(2, 2);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {1, 0};
  b.values = {Complex(1, 0), Complex(1, 0)};
  return MakeCase("no_explicit_zeros_in_output", a, b, [](const CRSMatrix &c) { return HasNoExplicitZeros(c); });
}(),

    // Non-trivial
    [] {
  CRSMatrix a(2, 3);
  a.row_ptr = {0, 2, 3};
  a.col_indices = {0, 2, 1};
  a.values = {Complex(1, 0), Complex(2, 1), Complex(3, 0)};

  CRSMatrix b(3, 4);
  b.row_ptr = {0, 1, 3, 4};
  b.col_indices = {1, 2, 3, 0};
  b.values = {Complex(1, 1), Complex(2, 0), Complex(1, 1), Complex(3, 0)};

  const std::vector<std::vector<Complex>> expected = {{Complex(6, 3), Complex(1, 1), Complex(0, 0), Complex(0, 0)},
                                                      {Complex(0, 0), Complex(0, 0), Complex(6, 0), Complex(3, 3)}};
  return MakeCase("rectangular_with_accumulation", a, b,
                  [expected](const CRSMatrix &c) { return DenseEqual(expected, c); });
}(), [] {
  CRSMatrix a(3, 2);
  a.row_ptr = {0, 2, 4, 6};
  a.col_indices = {0, 1, 0, 1, 0, 1};
  a.values = {Complex(1, 0), Complex(0, 1), Complex(2, 0), Complex(0, 2), Complex(3, 0), Complex(0, 3)};

  CRSMatrix b(2, 1);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 0};
  b.values = {Complex(1, 1), Complex(1, -1)};

  CRSMatrix expected(3, 1);
  expected.row_ptr = {0, 1, 2, 3};
  expected.col_indices = {0, 0, 0};
  expected.values = {Complex(2, 2), Complex(4, 4), Complex(6, 6)};
  return MakeCase("multiple_rows_contribute_to_same_column", a, b,
                  [expected](const CRSMatrix &c) { return CrsEqual(expected, c); });
}(), [] {
  CRSMatrix a(5, 5);
  a.row_ptr = {0, 1, 1, 1, 1, 2};
  a.col_indices = {0, 4};
  a.values = {Complex(1, 0), Complex(0, 1)};

  CRSMatrix b(5, 5);
  b.row_ptr = {0, 1, 1, 1, 1, 2};
  b.col_indices = {0, 4};
  b.values = {Complex(2, 0), Complex(0, -1)};
  return MakeCase("corner_elements_only_5x5", a, b, [](const CRSMatrix &c) { return CheckCornerElements5x5(c); });
}(), [] {
  CRSMatrix a(1, 4);
  a.row_ptr = {0, 4};
  a.col_indices = {0, 1, 2, 3};
  a.values = {Complex(1, 1), Complex(2, 0), Complex(3, -1), Complex(0, 4)};

  CRSMatrix i(4, 4);
  i.row_ptr = {0, 1, 2, 3, 4};
  i.col_indices = {0, 1, 2, 3};
  i.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0), Complex(1, 0)};
  return MakeCase("dense_row_times_identity", a, i, [](const CRSMatrix &c) { return CheckDenseRowTimesIdentity(c); });
}(), [] {
  CRSMatrix j(2, 2);
  j.row_ptr = {0, 1, 2};
  j.col_indices = {1, 0};
  j.values = {Complex(1, 0), Complex(-1, 0)};

  CRSMatrix minus_i(2, 2);
  minus_i.row_ptr = {0, 1, 2};
  minus_i.col_indices = {0, 1};
  minus_i.values = {Complex(-1, 0), Complex(-1, 0)};
  return MakeCase("matrix_squared_known_result", j, j, [minus_i](const CRSMatrix &c) { return CrsEqual(minus_i, c); });
}()};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<VidermanASparseMatrixMultCRSComplexSEQ, InType>(
                                               kTestParam, PPC_SETTINGS_viderman_a_sparse_matrix_mult_crs_complex),
                                           ppc::util::AddFuncTask<VidermanASparseMatrixMultCRSComplexOMP, InType>(
                                               kTestParam, PPC_SETTINGS_viderman_a_sparse_matrix_mult_crs_complex));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = VidermanRunFuncTests::PrintFuncTestName<VidermanRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(CRSComplexTests, VidermanRunFuncTests, kGtestValues, kTestName);

}  // namespace

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
