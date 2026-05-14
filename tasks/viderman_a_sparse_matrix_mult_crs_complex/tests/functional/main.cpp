#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/common/include/common.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/omp/include/ops_omp.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/seq/include/ops_seq.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/tbb/include/ops_tbb.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {
namespace {

constexpr double kTestTol = 1e-12;

// Структура параметров с принудительной инициализацией всех полей
struct MatrixTestParam {
  std::string name;
  CRSMatrix a;
  CRSMatrix b;
  std::function<bool(const CRSMatrix &)> check{nullptr};
  bool expect_valid = true;
};

void PrintTo(const ppc::util::FuncTestParam<InType, OutType, MatrixTestParam> &param, std::ostream *os) {
  *os << std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kNameTest)>(param);
}

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

// Хелперы для создания объектов параметров
MatrixTestParam MakeInvalid(const std::string &name, CRSMatrix a, CRSMatrix b) {
  return MatrixTestParam{.name = name, .a = std::move(a), .b = std::move(b), .check = nullptr, .expect_valid = false};
}

MatrixTestParam MakeCase(const std::string &name, CRSMatrix a, CRSMatrix b,
                         std::function<bool(const CRSMatrix &)> check) {
  return MatrixTestParam{
      .name = name, .a = std::move(a), .b = std::move(b), .check = std::move(check), .expect_valid = true};
}

testing::AssertionResult CheckInvalidCase(const ppc::util::FuncTestParam<InType, OutType, MatrixTestParam> &test_param,
                                          const InType &input_data) {
  auto task = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTaskGetter)>(test_param)(input_data);
  if (task->Validation()) {
    return testing::AssertionFailure() << "Validation() should return false";
  }
  try {
    task->Validation();
    return testing::AssertionFailure() << "Validation() should throw std::runtime_error";
  } catch (const std::runtime_error &) {
    return testing::AssertionSuccess();
  } catch (...) {
    return testing::AssertionFailure() << "Validation() threw an unexpected exception";
  }
}

}  // namespace

class VidermanRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, MatrixTestParam> {
 public:
  static std::string PrintTestParam(const testing::TestParamInfo<VidermanRunFuncTests::ParamType> &info) {
    return std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kNameTest)>(info.param) + "_" +
           std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(info.param).name;
  }

 protected:
  void SetUp() override {
    const auto &params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
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
  MatrixTestParam test_case_{};
};

namespace {

TEST_P(VidermanRunFuncTests, CRSComplexMult) {
  const auto &test_param = GetParam();
  const auto &test_case = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(test_param);

  if (!test_case.expect_valid) {
    EXPECT_TRUE(CheckInvalidCase(test_param, GetTestInputData()));
    return;
  }

  ExecuteTest(test_param);
}

const std::array<MatrixTestParam, 6> kTestParam = {
    MakeInvalid("incompatible_dimensions", CRSMatrix(2, 3), CRSMatrix(4, 5)),

    [] {
  CRSMatrix a(1, 1);
  a.row_ptr = {0, 1};
  a.col_indices = {0};
  a.values = {Complex(3.0, 4.0)};
  CRSMatrix b(1, 1);
  b.row_ptr = {0, 1};
  b.col_indices = {0};
  b.values = {Complex(1.0, -2.0)};
  CRSMatrix exp(1, 1);
  exp.row_ptr = {0, 1};
  exp.col_indices = {0};
  exp.values = {Complex(11.0, -2.0)};
  return MakeCase("single_element", std::move(a), std::move(b), [exp](const CRSMatrix &c) { return CrsEqual(exp, c); });
}(),

    MakeCase("both_zero_matrices", CRSMatrix(3, 4), CRSMatrix(4, 5),
             [](const CRSMatrix &c) { return c.rows == 3 && c.cols == 5 && c.values.empty() && c.IsValid(); }),

    [] {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 1};
  a.values = {Complex(1, 0), Complex(2, 0)};
  return MakeCase("nonzero_a_zero_b", std::move(a), CRSMatrix(2, 2),
                  [](const CRSMatrix &c) { return c.values.empty(); });
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
  return MakeCase("col_times_row", std::move(a), std::move(b),
                  [](const CRSMatrix &c) { return c.rows == 2 && c.cols == 2 && c.values.size() == 4; });
}(),

    MakeCase("identity_check", [] {
  CRSMatrix i(2, 2);
  i.row_ptr = {0, 1, 2};
  i.col_indices = {0, 1};
  i.values = {1, 1};
  return i;
}(), [] {
  CRSMatrix i(2, 2);
  i.row_ptr = {0, 1, 2};
  i.col_indices = {0, 1};
  i.values = {1, 1};
  return i;
}(), [](const CRSMatrix &c) { return c.values.size() == 2; })};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<VidermanASparseMatrixMultCRSComplexSEQ, InType>(
                                               kTestParam, PPC_SETTINGS_viderman_a_sparse_matrix_mult_crs_complex),
                                           ppc::util::AddFuncTask<VidermanASparseMatrixMultCRSComplexOMP, InType>(
                                               kTestParam, PPC_SETTINGS_viderman_a_sparse_matrix_mult_crs_complex),
                                           ppc::util::AddFuncTask<VidermanASparseMatrixMultCRSComplexTBB, InType>(
                                               kTestParam, PPC_SETTINGS_viderman_a_sparse_matrix_mult_crs_complex));

INSTANTIATE_TEST_SUITE_P(CRSComplexTests, VidermanRunFuncTests, ppc::util::ExpandToValues(kTestTasksList),
                         VidermanRunFuncTests::PrintTestParam);

}  // namespace
}  // namespace viderman_a_sparse_matrix_mult_crs_complex
