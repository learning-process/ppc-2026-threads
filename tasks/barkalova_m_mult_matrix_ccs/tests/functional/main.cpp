/*
#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// #include "barkalova_m_mult_matrix_ccs/all/include/ops_all.hpp"
#include "barkalova_m_mult_matrix_ccs/common/include/common.hpp"
// #include "barkalova_m_mult_matrix_ccs/omp/include/ops_omp.hpp"
#include "barkalova_m_mult_matrix_ccs/seq/include/ops_seq.hpp"
// #include "barkalova_m_mult_matrix_ccs/stl/include/ops_stl.hpp"
// #include "barkalova_m_mult_matrix_ccs/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace barkalova_m_mult_matrix_ccs {

class BarkalovaMMultMatrixCcsFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    int width = -1;
    int height = -1;
    int channels = -1;
    std::vector<uint8_t> img;
    // Read image
    {
      std::string abs_path = ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_barkalova_m_mult_matrix_ccs), "pic.jpg");
      auto *data = stbi_load(abs_path.c_str(), &width, &height, &channels, 0);
      if (data == nullptr) {
        throw std::runtime_error("Failed to load image: " + std::string(stbi_failure_reason()));
      }
      img = std::vector<uint8_t>(data, data + (static_cast<ptrdiff_t>(width * height * channels)));
      stbi_image_free(data);
      if (std::cmp_not_equal(width, height)) {
        throw std::runtime_error("width != height: ");
      }
    }

    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = width - height + std::min(std::accumulate(img.begin(), img.end(), 0), channels);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (input_data_ == output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = 0;
};

namespace {

TEST_P(BarkalovaMMultMatrixCcsFuncTestsThreads, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(3, "3"), std::make_tuple(5, "5"), std::make_tuple(7, "7")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<BarkalobaMMultMatrixCcsSEQ, InType>(kTestParam, PPC_SETTINGS_barkalova_m_mult_matrix_ccs));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    BarkalovaMMultMatrixCcsFuncTestsThreads::PrintFuncTestName<BarkalovaMMultMatrixCcsFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, BarkalovaMMultMatrixCcsFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace barkalova_m_mult_matrix_ccs
*/

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "barkalova_m_mult_matrix_ccs/common/include/common.hpp"
#include "barkalova_m_mult_matrix_ccs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace barkalova_m_mult_matrix_ccs {

class BarkalovaMatrixMultiplyFixedTest : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return "test_" + std::to_string(std::get<0>(test_param));
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int test_case = std::get<0>(params);

    CCSMatrix A, B;
    expected_result_ = CCSMatrix();

    switch (test_case) {
      case 1: {
        // Тест 1: Диагональные матрицы с комплексными числами
        A.rows = 2;
        A.cols = 2;
        A.col_ptrs = {0, 1, 2};
        A.row_indices = {0, 1};
        A.values = {Complex(1.0, 2.0), Complex(3.0, 4.0)};
        A.nnz = static_cast<int>(A.values.size());

        B.rows = 2;
        B.cols = 2;
        B.col_ptrs = {0, 1, 2};
        B.row_indices = {0, 1};
        B.values = {Complex(2.0, -1.0), Complex(1.0, 1.0)};
        B.nnz = static_cast<int>(B.values.size());

        expected_result_.rows = 2;
        expected_result_.cols = 2;
        expected_result_.col_ptrs = {0, 1, 2};
        expected_result_.row_indices = {0, 1};
        expected_result_.values = {Complex(4.0, 3.0), Complex(-1.0, 7.0)};
        expected_result_.nnz = static_cast<int>(expected_result_.values.size());
        break;
      }

      case 2: {
        // Тест 2: Единичная матрица с комплексными числами
        A.rows = 3;
        A.cols = 3;
        A.col_ptrs = {0, 1, 2, 3};
        A.row_indices = {0, 1, 2};
        A.values = {Complex(1.0, 0.0), Complex(1.0, 0.0), Complex(1.0, 0.0)};
        A.nnz = static_cast<int>(A.values.size());

        B.rows = 3;
        B.cols = 3;
        B.col_ptrs = {0, 1, 2, 3};
        B.row_indices = {0, 1, 2};
        B.values = {Complex(2.0, 3.0), Complex(4.0, -1.0), Complex(5.0, 2.0)};
        B.nnz = static_cast<int>(B.values.size());

        expected_result_.rows = 3;
        expected_result_.cols = 3;
        expected_result_.col_ptrs = {0, 1, 2, 3};
        expected_result_.row_indices = {0, 1, 2};
        expected_result_.values = {Complex(2.0, 3.0), Complex(4.0, -1.0), Complex(5.0, 2.0)};
        expected_result_.nnz = static_cast<int>(expected_result_.values.size());
        break;
      }

      case 3: {
        A.rows = 2;
        A.cols = 2;
        A.col_ptrs = {0, 1, 2};  // По одному элементу в каждом столбце
        A.row_indices = {0, 1};  // A[0][0] и A[1][1]
        A.values = {Complex(1.0, 0.0), Complex(2.0, 0.0)};
        A.nnz = static_cast<int>(A.values.size());

        B.rows = 2;
        B.cols = 2;
        B.col_ptrs = {0, 1, 2};  // По одному элементу в каждом столбце
        B.row_indices = {0, 1};  // B[0][0] и B[1][1]
        B.values = {Complex(3.0, 0.0), Complex(4.0, 0.0)};
        B.nnz = static_cast<int>(B.values.size());

        // Результат A × B:
        // C[0][0] = A[0][0]×B[0][0] = 1×3 = 3
        // C[1][1] = A[1][1]×B[1][1] = 2×4 = 8
        expected_result_.rows = 2;
        expected_result_.cols = 2;
        expected_result_.col_ptrs = {0, 1, 2};
        expected_result_.row_indices = {0, 1};
        expected_result_.values = {Complex(3.0, 0.0), Complex(8.0, 0.0)};
        expected_result_.nnz = 2;
        break;
      }

      case 4: {
        // Тест 4: Дополнительный тест с диагональными матрицами
        A.rows = 2;
        A.cols = 2;
        A.col_ptrs = {0, 1, 2};
        A.row_indices = {0, 1};
        A.values = {Complex(1.0, 1.0), Complex(2.0, -1.0)};
        A.nnz = static_cast<int>(A.values.size());

        B.rows = 2;
        B.cols = 2;
        B.col_ptrs = {0, 1, 2};
        B.row_indices = {0, 1};
        B.values = {Complex(3.0, 2.0), Complex(4.0, -3.0)};
        B.nnz = static_cast<int>(B.values.size());

        expected_result_.rows = 2;
        expected_result_.cols = 2;
        expected_result_.col_ptrs = {0, 1, 2};
        expected_result_.row_indices = {0, 1};
        expected_result_.values = {Complex(1.0, 5.0), Complex(5.0, -10.0)};
        expected_result_.nnz = static_cast<int>(expected_result_.values.size());
        break;
      }

      default:
        throw std::runtime_error("Unknown test case");
    }

    input_data_ = std::make_pair(A, B);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const double eps = 1e-10;

    if (output_data.rows != expected_result_.rows || output_data.cols != expected_result_.cols) {
      return false;
    }

    if (output_data.nnz != expected_result_.nnz) {
      return false;  // Убрали std::cout
    }

    if (output_data.col_ptrs.size() != expected_result_.col_ptrs.size()) {
      return false;
    }

    for (size_t i = 0; i < output_data.col_ptrs.size(); ++i) {
      if (output_data.col_ptrs[i] != expected_result_.col_ptrs[i]) {
        return false;  // Убрали std::cout
      }
    }

    if (output_data.row_indices.size() != expected_result_.row_indices.size()) {
      return false;
    }

    for (size_t i = 0; i < output_data.row_indices.size(); ++i) {
      if (output_data.row_indices[i] != expected_result_.row_indices[i]) {
        return false;  // Убрали std::cout
      }
    }

    for (size_t i = 0; i < output_data.values.size(); ++i) {
      if (std::abs(output_data.values[i].real() - expected_result_.values[i].real()) > eps ||
          std::abs(output_data.values[i].imag() - expected_result_.values[i].imag()) > eps) {
        return false;  // Убрали std::cout
      }
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  CCSMatrix expected_result_;
};

namespace {

TEST_P(BarkalovaMatrixMultiplyFixedTest, MatrixMultiplyFixedTest) {
  ExecuteTest(GetParam());
}

// Добавлен четвертый тест
const std::array<TestType, 4> kFixedTestParams = {std::make_tuple(1, ""), std::make_tuple(2, ""),
                                                  std::make_tuple(3, ""), std::make_tuple(4, "")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<BarkalovaMMultMatrixCcsSEQ, InType>(
    kFixedTestParams, PPC_SETTINGS_barkalova_m_mult_matrix_ccs));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = BarkalovaMatrixMultiplyFixedTest::PrintFuncTestName<BarkalovaMatrixMultiplyFixedTest>;

INSTANTIATE_TEST_SUITE_P(FixedMatrixTests, BarkalovaMatrixMultiplyFixedTest, kGtestValues, kTestName);

}  // namespace

}  // namespace barkalova_m_mult_matrix_ccs
