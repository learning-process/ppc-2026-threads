#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
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
#include "viderman_a_sparse_matrix_mult_crs_complex/tbb/include/ops_tbb.hpp"

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

CRSMatrix RunSeqMultiply(const CRSMatrix &a, const CRSMatrix &b) {
  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(a, b));
  if (!(task.Validation() && task.PreProcessing() && task.Run() && task.PostProcessing())) {
    return {};
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
    return test_case_.check ? test_case_.check(output_data) : false;
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

const std::array<TestCaseType, 15> kTestParam = {MakeInvalid("incompatible_dims", CRSMatrix(2, 3), CRSMatrix(4, 5)),
                                                 [] {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 10};
  a.values = {1, 2};
  return MakeInvalid("idx_out_of_range", a, CRSMatrix(2, 2));
}(),
                                                 [] {
  CRSMatrix a(1, 2);
  a.row_ptr = {0, 2};
  a.col_indices = {1, 0};
  a.values = {1, 2};
  return MakeInvalid("unsorted_indices", a, CRSMatrix(2, 2));
}(),

                                                 [] {
  CRSMatrix a(1, 1);
  a.row_ptr = {0, 1};
  a.col_indices = {0};
  a.values = {Complex(3, 4)};
  CRSMatrix b(1, 1);
  b.row_ptr = {0, 1};
  b.col_indices = {0};
  b.values = {Complex(1, -2)};
  CRSMatrix exp(1, 1);
  exp.row_ptr = {0, 1};
  exp.col_indices = {0};
  exp.values = {Complex(11, -2)};
  return MakeCase("single_element", a, b, [exp](const CRSMatrix &c) { return CrsEqual(exp, c); });
}(), MakeCase("both_empty", CRSMatrix(2, 2), CRSMatrix(2, 2), [](const CRSMatrix &c) { return c.values.empty(); }),
                                                 [] {
  CRSMatrix a(1, 2);
  a.row_ptr = {0, 2};
  a.col_indices = {0, 1};
  a.values = {1, -1};
  CRSMatrix b(2, 1);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 0};
  b.values = {Complex(0, 1), Complex(0, 1)};
  return MakeCase("cancellation_to_zero", a, b, [](const CRSMatrix &c) { return c.values.empty(); });
}(),
                                                 [] {
  CRSMatrix a(2, 1);
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 0};
  a.values = {Complex(1, 1), 2};
  CRSMatrix b(1, 2);
  b.row_ptr = {0, 2};
  b.col_indices = {0, 1};
  b.values = {3, Complex(0, 1)};
  return MakeCase("col_times_row", a, b, [](const CRSMatrix &c) { return DenseEqual(kDense2x2Expected, c); });
}(),

                                                 // --- Functional (8) ---
                                                 [] {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 1};
  a.values = {Complex(0, 1), Complex(0, 1)};
  CRSMatrix exp(2, 2);
  exp.row_ptr = {0, 1, 2};
  exp.col_indices = {0, 1};
  exp.values = {-1, -1};
  return MakeCase("imaginary_square", a, a, [exp](const CRSMatrix &c) { return CrsEqual(exp, c); });
}(), [] {
  CRSMatrix a(3, 3);
  a.row_ptr = {0, 2, 3, 4};
  a.col_indices = {0, 2, 1, 0};
  a.values = {Complex(1, 2), 3, Complex(0, 1), Complex(5, -1)};
  CRSMatrix i(3, 3);
  i.row_ptr = {0, 1, 2, 3};
  i.col_indices = {0, 1, 2};
  i.values = {1, 1, 1};
  return MakeCase("identity_right", a, i, [a](const CRSMatrix &c) { return CrsEqual(a, c); });
}(), [] {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 2, 3};
  a.col_indices = {0, 1, 0};
  a.values = {Complex(1, 1), 2, Complex(0, 1)};
  CRSMatrix b(2, 2);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 1};
  b.values = {3, Complex(0, 2)};
  CRSMatrix c(2, 2);
  c.row_ptr = {0, 2, 2};
  c.col_indices = {0, 1};
  c.values = {1, Complex(1, 1)};
  return MakeCase("associativity", a, b, [a, b, c](const CRSMatrix &ab) {
    return DenseEqual(ToDense(RunSeqMultiply(ab, c)), RunSeqMultiply(a, RunSeqMultiply(b, c)));
  });
}(), [] {
  CRSMatrix p(2, 2);
  p.row_ptr = {0, 1, 2};
  p.col_indices = {1, 0};
  p.values = {1, 1};
  return MakeCase("permutation_square", p, p,
                  [](const CRSMatrix &c) { return c.values.size() == 2 && c.col_indices[0] == 0; });
}(), [] {
  CRSMatrix a(2, 3);
  a.row_ptr = {0, 2, 3};
  a.col_indices = {0, 2, 1};
  a.values = {1, Complex(2, 1), 3};
  CRSMatrix b(3, 2);
  b.row_ptr = {0, 1, 2, 3};
  b.col_indices = {0, 1, 0};
  b.values = {1, 1, 1};
  return MakeCase("rect_mult", a, b, [](const CRSMatrix &c) { return c.rows == 2 && c.cols == 2; });
}(), [] {
  CRSMatrix a(3, 3);
  a.row_ptr = {0, 1, 2, 3};
  a.col_indices = {0, 1, 2};
  a.values = {Complex(1, 1), 2, 3};
  CRSMatrix b(3, 3);
  b.row_ptr = {0, 1, 2, 3};
  b.col_indices = {0, 1, 2};
  b.values = {2, 2, 2};
  return MakeCase("diagonal", a, b,
                  [](const CRSMatrix &c) { return c.values.size() == 3 && ComplexNear(c.values[0], Complex(2, 2)); });
}(), [] {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 2, 4};
  a.col_indices = {0, 1, 0, 1};
  a.values = {1, 1, 1, 1};
  CRSMatrix b(2, 2);
  b.row_ptr = {0, 2, 4};
  b.col_indices = {0, 1, 0, 1};
  b.values = {1, 1, 1, 1};
  return MakeCase("dense_mult", a, b, [](const CRSMatrix &c) { return ComplexNear(c.values[0], 2.0); });
}(), [] {
  CRSMatrix a(1, 1);
  a.row_ptr = {0, 1};
  a.col_indices = {0};
  a.values = {Complex(3, 4)};
  CRSMatrix b(1, 1);
  b.row_ptr = {0, 1};
  b.col_indices = {0};
  b.values = {Complex(3, -4)};
  return MakeCase("conjugate", a, b, [](const CRSMatrix &c) { return ComplexNear(c.values[0], 25.0); });
}()};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<VidermanASparseMatrixMultCRSComplexSEQ, InType>(
                                               kTestParam, PPC_SETTINGS_viderman_a_sparse_matrix_mult_crs_complex),
                                           ppc::util::AddFuncTask<VidermanASparseMatrixMultCRSComplexOMP, InType>(
                                               kTestParam, PPC_SETTINGS_viderman_a_sparse_matrix_mult_crs_complex),
                                           ppc::util::AddFuncTask<VidermanASparseMatrixMultCRSComplexTBB, InType>(
                                               kTestParam, PPC_SETTINGS_viderman_a_sparse_matrix_mult_crs_complex));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kTestName = VidermanRunFuncTests::PrintFuncTestName<VidermanRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(CRSComplexTests, VidermanRunFuncTests, kGtestValues, kTestName);

}  // namespace
}  // namespace viderman_a_sparse_matrix_mult_crs_complex
