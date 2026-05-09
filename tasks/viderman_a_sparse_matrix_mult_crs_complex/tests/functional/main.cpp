#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <tuple>

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

  for (size_t i = 0; i < expected.values.size(); ++i) {
    if (!ComplexNear(expected.values[i], actual.values[i], tol)) {
      return false;
    }
  }

  return true;
}

struct MatrixTestParam {
  std::string name;
  CRSMatrix a;
  CRSMatrix b;

  std::function<bool(const CRSMatrix &)> check;
};

void ValidateInput(BaseTask &test_task) {
  ASSERT_TRUE(test_task.Validation());
}

void PrepareInput(BaseTask &test_task) {
  ASSERT_TRUE(test_task.PreProcessing());
}

void RunTestsuit(BaseTask &test_task) {
  ASSERT_TRUE(test_task.Run());
  ASSERT_TRUE(test_task.PostProcessing());
}

void CheckOutput(const MatrixTestParam &param, BaseTask &test_task) {
  OutType output = test_task.GetOutput();
  ASSERT_TRUE(param.check(output));
}

void RunMatrixMultTest(const MatrixTestParam &param, BaseTask &test_task) {
  ValidateInput(test_task);
  PrepareInput(test_task);
  RunTestsuit(test_task);
  CheckOutput(param, test_task);
}

}  // namespace

class VidermanRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, MatrixTestParam> {
 public:
  static std::string PrintTestParam(const MatrixTestParam &test_param) {
    return test_param.name;
  }

 protected:
  void SetUp() override {
    test_param_ = std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return test_param_.check ? test_param_.check(output_data) : false;
  }

  InType GetTestInputData() final {
    return std::make_tuple(test_param_.a, test_param_.b);
  }

 private:
  MatrixTestParam test_param_;
};

namespace {

TEST_P(VidermanRunFuncTests, CRSComplexMultFuncTest) {
  const auto &param = std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
  auto get_task = std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTaskGetter)>(GetParam());

  auto ptr_task = get_task(std::make_tuple(param.a, param.b));

  RunMatrixMultTest(param, *ptr_task);
}

const std::array<MatrixTestParam, 12> kTestParam = {
    [] {
  CRSMatrix a(1, 1);
  CRSMatrix b(1, 1);
  CRSMatrix exp(1, 1);

  a.row_ptr = {0, 1};
  a.col_indices = {0};
  a.values = {Complex(3, 4)};

  b.row_ptr = {0, 1};
  b.col_indices = {0};
  b.values = {Complex(1, -2)};

  exp.row_ptr = {0, 1};
  exp.col_indices = {0};
  exp.values = {Complex(11, -2)};

  return MatrixTestParam{.name = "single_element",
                         .a = a,
                         .b = b,

                         .check = [exp](const CRSMatrix &c) { return CrsEqual(exp, c); }};
}(),

    [] {
  CRSMatrix a(1, 2);
  CRSMatrix b(2, 1);

  a.row_ptr = {0, 2};
  a.col_indices = {0, 1};
  a.values = {1, -1};

  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 0};
  b.values = {Complex(0, 1), Complex(0, 1)};

  return MatrixTestParam{.name = "cancellation_to_zero",
                         .a = a,
                         .b = b,

                         .check = [](const CRSMatrix &c) { return c.values.empty(); }};
}(),

    [] {
  CRSMatrix a(2, 2);
  CRSMatrix exp(2, 2);

  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 1};
  a.values = {Complex(0, 1), Complex(0, 1)};

  exp.row_ptr = {0, 1, 2};
  exp.col_indices = {0, 1};
  exp.values = {-1, -1};

  return MatrixTestParam{.name = "imaginary_square",
                         .a = a,
                         .b = a,

                         .check = [exp](const CRSMatrix &c) { return CrsEqual(exp, c); }};
}(),

    [] {
  CRSMatrix a(3, 3);
  CRSMatrix i(3, 3);

  a.row_ptr = {0, 2, 3, 4};
  a.col_indices = {0, 2, 1, 0};
  a.values = {Complex(1, 2), 3, Complex(0, 1), Complex(5, -1)};

  i.row_ptr = {0, 1, 2, 3};
  i.col_indices = {0, 1, 2};
  i.values = {1, 1, 1};

  return MatrixTestParam{.name = "identity_right",
                         .a = a,
                         .b = i,

                         .check = [a](const CRSMatrix &c) { return CrsEqual(a, c); }};
}(),

    [] {
  CRSMatrix p(2, 2);

  p.row_ptr = {0, 1, 2};
  p.col_indices = {1, 0};
  p.values = {1, 1};

  return MatrixTestParam{.name = "permutation_square",
                         .a = p,
                         .b = p,

                         .check = [](const CRSMatrix &c) { return c.values.size() == 2 && c.col_indices[0] == 0; }};
}(),

    [] {
  CRSMatrix a(2, 3);
  CRSMatrix b(3, 2);

  a.row_ptr = {0, 2, 3};
  a.col_indices = {0, 2, 1};
  a.values = {1, Complex(2, 1), 3};

  b.row_ptr = {0, 1, 2, 3};
  b.col_indices = {0, 1, 0};
  b.values = {1, 1, 1};

  return MatrixTestParam{.name = "rect_mult",
                         .a = a,
                         .b = b,

                         .check = [](const CRSMatrix &c) { return c.rows == 2 && c.cols == 2; }};
}(),

    [] {
  CRSMatrix a(3, 3);
  CRSMatrix b(3, 3);

  a.row_ptr = {0, 1, 2, 3};
  a.col_indices = {0, 1, 2};
  a.values = {Complex(1, 1), 2, 3};

  b.row_ptr = {0, 1, 2, 3};
  b.col_indices = {0, 1, 2};
  b.values = {2, 2, 2};

  return MatrixTestParam{
      .name = "diagonal",
      .a = a,
      .b = b,

      .check = [](const CRSMatrix &c) { return c.values.size() == 3 && ComplexNear(c.values[0], Complex(2, 2)); }};
}(),

    [] {
  CRSMatrix a(2, 2);
  CRSMatrix b(2, 2);

  a.row_ptr = {0, 2, 4};
  a.col_indices = {0, 1, 0, 1};
  a.values = {1, 1, 1, 1};

  b.row_ptr = {0, 2, 4};
  b.col_indices = {0, 1, 0, 1};
  b.values = {1, 1, 1, 1};

  return MatrixTestParam{.name = "dense_mult",
                         .a = a,
                         .b = b,

                         .check = [](const CRSMatrix &c) { return ComplexNear(c.values[0], 2.0); }};
}(),

    [] {
  CRSMatrix a(1, 1);
  CRSMatrix b(1, 1);

  a.row_ptr = {0, 1};
  a.col_indices = {0};
  a.values = {Complex(3, 4)};

  b.row_ptr = {0, 1};
  b.col_indices = {0};
  b.values = {Complex(3, -4)};

  return MatrixTestParam{.name = "conjugate",
                         .a = a,
                         .b = b,

                         .check = [](const CRSMatrix &c) { return ComplexNear(c.values[0], 25.0); }};
}(),

    MatrixTestParam{.name = "both_empty",
                    .a = CRSMatrix(2, 2),
                    .b = CRSMatrix(2, 2),

                    .check = [](const CRSMatrix &c) { return c.values.empty(); }},

    MatrixTestParam{.name = "zero_1x1",
                    .a = CRSMatrix(1, 1),
                    .b = CRSMatrix(1, 1),

                    .check = [](const CRSMatrix &c) { return c.values.empty(); }},

    MatrixTestParam{.name = "sparse_large_empty",
                    .a = CRSMatrix(100, 100),
                    .b = CRSMatrix(100, 100),

                    .check = [](const CRSMatrix &c) { return c.values.empty(); }}};

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
