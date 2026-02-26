#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "kurpiakov_a_sp_comp_mat_mul/common/include/common.hpp"
#include "kurpiakov_a_sp_comp_mat_mul/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kurpiakov_a_sp_comp_mat_mul {

class KurpiakovRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType& test_param) {
    return std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_output_ = std::get<2>(params);
  }

  bool CheckTestOutputData(OutType& output_data) final {
    return expected_output_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_output_;
};

namespace {

// Test 1: Identity * A = A
InType MakeIdentityTest() {
  SparseMatrix eye(2, 2, {ComplexD(1, 0), ComplexD(1, 0)}, {0, 1}, {0, 1, 2});
  SparseMatrix a(2, 2, {ComplexD(1, 0), ComplexD(2, 1), ComplexD(3, 0)}, {0, 1, 1}, {0, 2, 3});
  return {eye, a};
}
OutType MakeIdentityExpected() {
  return {2, 2, {ComplexD(1, 0), ComplexD(2, 1), ComplexD(3, 0)}, {0, 1, 1}, {0, 2, 3}};
}

// Test 2: A * Zero = Zero
InType MakeZeroTest() {
  SparseMatrix a(2, 2, {ComplexD(1, 0), ComplexD(2, 1)}, {0, 1}, {0, 1, 2});
  SparseMatrix zero(2, 2);
  return {a, zero};
}
OutType MakeZeroExpected() {
  return {2, 2};
}

// Test 3: Diagonal * Diagonal
InType MakeDiagTest() {
  SparseMatrix a(2, 2, {ComplexD(1, 1), ComplexD(2, 0)}, {0, 1}, {0, 1, 2});
  SparseMatrix b(2, 2, {ComplexD(2, 0), ComplexD(0, 1)}, {0, 1}, {0, 1, 2});
  return {a, b};
}
OutType MakeDiagExpected() {
  return {2, 2, {ComplexD(2, 2), ComplexD(0, 2)}, {0, 1}, {0, 1, 2}};
}

// Test 4: General 2x2 complex
InType MakeGeneralTest() {
  SparseMatrix a(2, 2, {ComplexD(1, 1), ComplexD(1, -1)}, {0, 1}, {0, 1, 2});
  SparseMatrix b(2, 2, {ComplexD(2, 0), ComplexD(3, 1)}, {0, 1}, {0, 1, 2});
  return {a, b};
}
OutType MakeGeneralExpected() {
  return {2, 2, {ComplexD(2, 2), ComplexD(4, -2)}, {0, 1}, {0, 1, 2}};
}

// Test 5: Scalar i * i = -1
InType MakeScalarITest() {
  SparseMatrix si(1, 1, {ComplexD(0, 1)}, {0}, {0, 1});
  return {si, si};
}
OutType MakeScalarIExpected() {
  return {1, 1, {ComplexD(-1, 0)}, {0}, {0, 1}};
}

// Test 6: Scalar (3+4i)*(1-2i) = 11-2i
InType MakeScalarComplexTest() {
  SparseMatrix a(1, 1, {ComplexD(3, 4)}, {0}, {0, 1});
  SparseMatrix b(1, 1, {ComplexD(1, -2)}, {0}, {0, 1});
  return {a, b};
}
OutType MakeScalarComplexExpected() {
  return {1, 1, {ComplexD(11, -2)}, {0}, {0, 1}};
}

// Test 7: Rectangular A(2x3) * B(3x2)
InType MakeRectTest() {
  SparseMatrix a(2, 3, {ComplexD(1, 0), ComplexD(2, 0), ComplexD(1, 0)}, {0, 2, 1}, {0, 2, 3});
  SparseMatrix b(3, 2, {ComplexD(1, 0), ComplexD(1, 0), ComplexD(1, 0), ComplexD(1, 0)}, {0, 1, 0, 1}, {0, 2, 3, 4});
  return {a, b};
}
OutType MakeRectExpected() {
  return {2, 2, {ComplexD(1, 0), ComplexD(3, 0), ComplexD(1, 0)}, {0, 1, 0}, {0, 2, 3}};
}

// Test 8: Both empty 3x3
InType MakeEmptyTest() {
  SparseMatrix empty(3, 3);
  return {empty, empty};
}
OutType MakeEmptyExpected() {
  return {3, 3};
}

// Test 9: Dense 2x2
InType MakeDenseTest() {
  SparseMatrix a(2, 2, {ComplexD(1, 0), ComplexD(2, 0), ComplexD(3, 0), ComplexD(4, 0)}, {0, 1, 0, 1}, {0, 2, 4});
  SparseMatrix b(2, 2, {ComplexD(5, 0), ComplexD(6, 0), ComplexD(7, 0), ComplexD(8, 0)}, {0, 1, 0, 1}, {0, 2, 4});
  return {a, b};
}
OutType MakeDenseExpected() {
  return {2, 2, {ComplexD(19, 0), ComplexD(22, 0), ComplexD(43, 0), ComplexD(50, 0)}, {0, 1, 0, 1}, {0, 2, 4}};
}

// Test 10: Row * Column = scalar
InType MakeRowColTest() {
  SparseMatrix a(1, 3, {ComplexD(1, 0), ComplexD(2, 0), ComplexD(3, 0)}, {0, 1, 2}, {0, 3});
  SparseMatrix b(3, 1, {ComplexD(4, 0), ComplexD(5, 0), ComplexD(6, 0)}, {0, 0, 0}, {0, 1, 2, 3});
  return {a, b};
}
OutType MakeRowColExpected() {
  return {1, 1, {ComplexD(32, 0)}, {0}, {0, 1}};
}

// Test 11: Column * Row = outer product 3x3
InType MakeOuterProductTest() {
  SparseMatrix a(3, 1, {ComplexD(1, 0), ComplexD(2, 0), ComplexD(3, 0)}, {0, 0, 0}, {0, 1, 2, 3});
  SparseMatrix b(1, 3, {ComplexD(4, 0), ComplexD(5, 0), ComplexD(6, 0)}, {0, 1, 2}, {0, 3});
  return {a, b};
}
OutType MakeOuterProductExpected() {
  return {3,
          3,
          {ComplexD(4, 0), ComplexD(5, 0), ComplexD(6, 0), ComplexD(8, 0), ComplexD(10, 0), ComplexD(12, 0),
           ComplexD(12, 0), ComplexD(15, 0), ComplexD(18, 0)},
          {0, 1, 2, 0, 1, 2, 0, 1, 2},
          {0, 3, 6, 9}};
}

// Test 12: Sparse 3x3 with few nonzeros
InType MakeSparse3x3Test() {
  SparseMatrix a(3, 3, {ComplexD(1, 0), ComplexD(2, 0)}, {2, 0}, {0, 1, 1, 2});
  SparseMatrix b(3, 3, {ComplexD(1, 0), ComplexD(3, 0)}, {0, 2}, {0, 1, 1, 2});
  return {a, b};
}
OutType MakeSparse3x3Expected() {
  return {3, 3, {ComplexD(3, 0), ComplexD(2, 0)}, {2, 0}, {0, 1, 1, 2}};
}

// Test 13: Complex dense 2x2
InType MakeComplexDenseTest() {
  SparseMatrix a(2, 2, {ComplexD(1, 2), ComplexD(3, 4), ComplexD(9, 10), ComplexD(11, 12)}, {0, 1, 0, 1}, {0, 2, 4});
  SparseMatrix b(2, 2, {ComplexD(5, 6), ComplexD(7, 8), ComplexD(13, 14), ComplexD(15, 16)}, {0, 1, 0, 1}, {0, 2, 4});
  return {a, b};
}
OutType MakeComplexDenseExpected() {
  return {
      2, 2, {ComplexD(-24, 110), ComplexD(-28, 130), ComplexD(-40, 414), ComplexD(-44, 498)}, {0, 1, 0, 1}, {0, 2, 4}};
}

// Test 14: A * Identity = A (right multiply)
InType MakeRightIdentityTest() {
  SparseMatrix a(2, 2, {ComplexD(5, -3), ComplexD(7, 2), ComplexD(1, 1)}, {0, 0, 1}, {0, 1, 3});
  SparseMatrix eye(2, 2, {ComplexD(1, 0), ComplexD(1, 0)}, {0, 1}, {0, 1, 2});
  return {a, eye};
}
OutType MakeRightIdentityExpected() {
  return {2, 2, {ComplexD(5, -3), ComplexD(7, 2), ComplexD(1, 1)}, {0, 0, 1}, {0, 1, 3}};
}

// Test 15: Tridiagonal 4x4 * Identity
InType MakeTridiag4Test() {
  SparseMatrix a(4, 4,
                 {ComplexD(2, 0), ComplexD(1, 0), ComplexD(1, 0), ComplexD(2, 0), ComplexD(1, 0), ComplexD(1, 0),
                  ComplexD(2, 0), ComplexD(1, 0), ComplexD(1, 0), ComplexD(2, 0)},
                 {0, 1, 0, 1, 2, 1, 2, 3, 2, 3}, {0, 2, 5, 8, 10});
  SparseMatrix eye(4, 4, {ComplexD(1, 0), ComplexD(1, 0), ComplexD(1, 0), ComplexD(1, 0)}, {0, 1, 2, 3},
                   {0, 1, 2, 3, 4});
  return {a, eye};
}
OutType MakeTridiag4Expected() {
  return {4,
          4,
          {ComplexD(2, 0), ComplexD(1, 0), ComplexD(1, 0), ComplexD(2, 0), ComplexD(1, 0), ComplexD(1, 0),
           ComplexD(2, 0), ComplexD(1, 0), ComplexD(1, 0), ComplexD(2, 0)},
          {0, 1, 0, 1, 2, 1, 2, 3, 2, 3},
          {0, 2, 5, 8, 10}};
}

const std::array<TestType, 15> kTestParam = {
    std::make_tuple(MakeIdentityTest(), "identity_2x2", MakeIdentityExpected()),
    std::make_tuple(MakeZeroTest(), "zero_2x2", MakeZeroExpected()),
    std::make_tuple(MakeDiagTest(), "diag_2x2", MakeDiagExpected()),
    std::make_tuple(MakeGeneralTest(), "general_2x2", MakeGeneralExpected()),
    std::make_tuple(MakeScalarITest(), "scalar_i_squared", MakeScalarIExpected()),
    std::make_tuple(MakeScalarComplexTest(), "scalar_complex", MakeScalarComplexExpected()),
    std::make_tuple(MakeRectTest(), "rect_2x3_times_3x2", MakeRectExpected()),
    std::make_tuple(MakeEmptyTest(), "both_empty_3x3", MakeEmptyExpected()),
    std::make_tuple(MakeDenseTest(), "dense_2x2", MakeDenseExpected()),
    std::make_tuple(MakeRowColTest(), "row_times_col", MakeRowColExpected()),
    std::make_tuple(MakeOuterProductTest(), "outer_product", MakeOuterProductExpected()),
    std::make_tuple(MakeSparse3x3Test(), "sparse_3x3", MakeSparse3x3Expected()),
    std::make_tuple(MakeComplexDenseTest(), "complex_dense_2x2", MakeComplexDenseExpected()),
    std::make_tuple(MakeRightIdentityTest(), "right_identity_2x2", MakeRightIdentityExpected()),
    std::make_tuple(MakeTridiag4Test(), "tridiag_4x4_times_eye", MakeTridiag4Expected()),
};

TEST_P(KurpiakovRunFuncTestsThreads, SparseMatMulFromParams) {
  ExecuteTest(GetParam());
}

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<KurpiskovACRSMatMulSEQ, InType>(kTestParam, PPC_SETTINGS_kurpiakov_a_sp_comp_mat_mul));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kPerfTestName = KurpiakovRunFuncTestsThreads::PrintFuncTestName<KurpiakovRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(SparseMatMulSeqTests, KurpiakovRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kurpiakov_a_sp_comp_mat_mul
