#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "posternak_a_crs_mul_complex_matrix/common/include/common.hpp"
#include "posternak_a_crs_mul_complex_matrix/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace posternak_a_crs_mul_complex_matrix {

// CRS матрица из плотной(обычной)
CRSMatrix MakeCRS(int rows, int cols, const std::vector<std::vector<std::complex<double>>> &default_matrix) {
  CRSMatrix m;
  m.rows = rows;
  m.cols = cols;
  m.index_row.push_back(0);

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (std::abs(default_matrix[i][j]) > 1e-12) {
        m.values.push_back(default_matrix[i][j]);
        m.index_col.push_back(j);
      }
    }
    m.index_row.push_back(static_cast<int>(m.values.size()));
  }
  return m;
}

// сравнение с учетом погрешности
bool MatricesEqual(const CRSMatrix &a, const CRSMatrix &b) {
  if (a.rows != b.rows || a.cols != b.cols) {
    return false;
  }
  if (a.values.size() != b.values.size()) {
    return false;
  }

  for (size_t i = 0; i < a.values.size(); ++i) {
    if (std::abs(a.values[i] - b.values[i]) > 1e-9) {
      return false;
    }
    if (a.index_col[i] != b.index_col[i]) {
      return false;
    }
  }
  for (size_t i = 0; i < a.index_row.size(); ++i) {
    if (a.index_row[i] != b.index_row[i]) {
      return false;
    }
  }
  return true;
}

class PosternakARunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int test_id = std::get<0>(params);

    CRSMatrix &a = std::get<0>(input_data_);
    CRSMatrix &b = std::get<1>(input_data_);
    CRSMatrix &c = expected_result_;

    // A(2×3) × B(3×2) = C(2×2)
    if (test_id == 0) {
      std::vector<std::vector<std::complex<double>>> mat_a = {{{1, 0}, {0, 0}, {2, 0}}, {{0, 0}, {3, 0}, {0, 0}}};

      std::vector<std::vector<std::complex<double>>> mat_b = {{{4, 0}, {0, 0}}, {{0, 0}, {6, 0}}, {{5, 0}, {0, 0}}};

      std::vector<std::vector<std::complex<double>>> mat_c = {{{14, 0}, {0, 0}}, {{0, 0}, {18, 0}}};

      a = MakeCRS(2, 3, mat_a);
      b = MakeCRS(3, 2, mat_b);
      c = MakeCRS(2, 2, mat_c);
    }

    if (test_id == 1) {
      std::vector<std::vector<std::complex<double>>> mat_a = {{{1, 1}}};   // [1+i]
      std::vector<std::vector<std::complex<double>>> mat_b = {{{1, -1}}};  // [1-i]
      std::vector<std::vector<std::complex<double>>> mat_c = {{{2, 0}}};   // [2+0i]

      a = MakeCRS(1, 1, mat_a);
      b = MakeCRS(1, 1, mat_b);
      c = MakeCRS(1, 1, mat_c);
    }

    if (test_id == 2) {
      std::vector<std::vector<std::complex<double>>> mat_a = {
          {{1, 2}, {0, 0}, {0, 0}}, {{0, 0}, {0, 0}, {3, -1}}, {{0, 0}, {0, 1}, {0, 0}}};

      std::vector<std::vector<std::complex<double>>> mat_ones = {
          {{1, 0}, {0, 0}, {0, 0}}, {{0, 0}, {1, 0}, {0, 0}}, {{0, 0}, {0, 0}, {1, 0}}};

      a = MakeCRS(3, 3, mat_a);
      b = MakeCRS(3, 3, mat_ones);
      c = MakeCRS(3, 3, mat_a);
    }

    if (test_id == 3) {
      std::vector<std::vector<std::complex<double>>> empty_2x2(2, std::vector<std::complex<double>>(2, {0, 0}));

      a = MakeCRS(2, 2, empty_2x2);
      b = MakeCRS(2, 2, empty_2x2);
      c = MakeCRS(2, 2, empty_2x2);
    }

    if (test_id == 4) {
      std::vector<std::vector<std::complex<double>>> mat_a = {{{5, -3}, {0, 0}}, {{0, 0}, {1, 1}}};
      std::vector<std::vector<std::complex<double>>> zero_matrix(2, std::vector<std::complex<double>>(2, {0, 0}));

      a = MakeCRS(2, 2, mat_a);
      b = MakeCRS(2, 2, zero_matrix);
      c = MakeCRS(2, 2, zero_matrix);
    }

    if (test_id == 5) {
      // A(4×3) × B(3×4) = C(4×4)
      std::vector<std::vector<std::complex<double>>> mat_a = {
          {{1, 1}, {0, 0}, {2, -1}}, {{0, 0}, {3, 2}, {0, 0}}, {{0, 0}, {0, 0}, {1, 1}}, {{2, 0}, {0, 0}, {0, 0}}};

      std::vector<std::vector<std::complex<double>>> mat_b = {
          {{1, 0}, {0, 0}, {0, 1}, {0, 0}}, {{0, 0}, {2, 0}, {0, 0}, {1, -1}}, {{1, -1}, {0, 0}, {2, 0}, {0, 0}}};

      // C[0][0] = (1+i)*1 + (2-i)*(1-i) = 2-2i
      // C[0][2] = (1+i)*i + (2-i)*2 = 3-i
      // C[1][1] = (3+2i)*2 = 6+4i
      // C[1][3] = (3+2i)*(1-i) = 5-i
      // C[2][0] = (1+i)*(1-i) = 2
      // C[2][2] = (1+i)*2 = 2+2i
      // C[3][0] = 2*1 = 2
      // C[3][2] = 2*i = 2i
      std::vector<std::vector<std::complex<double>>> mat_c = {{{2, -2}, {0, 0}, {3, -1}, {0, 0}},
                                                              {{0, 0}, {6, 4}, {0, 0}, {5, -1}},
                                                              {{2, 0}, {0, 0}, {2, 2}, {0, 0}},
                                                              {{2, 0}, {0, 0}, {0, 2}, {0, 0}}};

      a = MakeCRS(4, 3, mat_a);
      b = MakeCRS(3, 4, mat_b);
      c = MakeCRS(4, 4, mat_c);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return MatricesEqual(expected_result_, output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_result_;
};

namespace {

TEST_P(PosternakARunFuncTestsThreads, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 6> kTestParams = {
    std::make_tuple(0, "default_test"),         std::make_tuple(1, "complex_1x1"),
    std::make_tuple(2, "multiply_by_identity"), std::make_tuple(3, "empty_matrices"),
    std::make_tuple(4, "multiply_by_zero"),     std::make_tuple(5, "complex_4x3_3x4")};

const auto kTestTasksList = ppc::util::AddFuncTask<PosternakACRSMulComplexMatrixSEQ, InType>(
    kTestParams, PPC_SETTINGS_posternak_a_crs_mul_complex_matrix);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = PosternakARunFuncTestsThreads::PrintFuncTestName<PosternakARunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(FuncTests, PosternakARunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace posternak_a_crs_mul_complex_matrix
