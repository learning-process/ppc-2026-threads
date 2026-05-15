#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <tuple>

#include "morozova_s_strassen_multiplication/common/include/common.hpp"
#include "morozova_s_strassen_multiplication/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace morozova_s_strassen_multiplication {

class MorozovaSStrassenMultiplicationFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    test_number_ = std::get<0>(params);
    SetupTestData();
  }

  bool CheckTestOutputData(OutType &output_data) override {
    if (test_number_ == 6 || test_number_ == 7) {
      return output_data.empty();
    }
    return ValidateMultiplicationResult(output_data);
  }

  InType GetTestInputData() override {
    return input_data_;
  }

 private:
  void SetupTestData() {
    switch (test_number_) {
      case 1:
        SetupTest1();
        break;
      case 2:
        SetupTest2();
        break;
      case 3:
        SetupTest3();
        break;
      case 4:
        SetupTest4();
        break;
      case 5:
        SetupTest5();
        break;
      case 6:
        SetupTest6();
        break;
      case 7:
        SetupTest7();
        break;
      default:
        SetupDefaultTest();
        break;
    }
  }

  void SetupTest1() {
    input_data_ = {2.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
  }

  void SetupTest2() {
    input_data_ = {4.0};
    AddIdentityMatrix(4);
    AddSequentialMatrix(4);
  }

  void SetupTest3() {
    input_data_ = {8.0};
    AddWeightedMatrix(8, [](int i, int j) { return static_cast<double>(i + 1) * (j + 1) * 0.5; });
    AddWeightedMatrix(8, [](int i, int j) { return static_cast<double>(i + j + 1) * 0.3; });
  }

  void SetupTest4() {
    input_data_ = {16.0};
    AddWeightedMatrix(16, [](int i, int j) { return std::sin(static_cast<double>(i + j)) * 10.0; });
    AddWeightedMatrix(16, [](int i, int j) { return std::cos(static_cast<double>(i - j)) * 5.0; });
  }

  void SetupTest5() {
    input_data_ = {32.0};
    AddWeightedMatrix(32, [](int i, int j) { return static_cast<double>((i * 32) + j + 1); });
    AddWeightedMatrix(32, [](int i, int j) { return static_cast<double>(((i + j) * 2) + 1); });
  }

  void SetupTest6() {
    input_data_ = {};
  }

  void SetupTest7() {
    input_data_ = {0.0, 1.0, 2.0, 3.0, 4.0};
  }

  void SetupDefaultTest() {
    input_data_ = {2.0, 1.0, 0.0, 0.0, 1.0, 1.0, 2.0, 3.0, 4.0};
  }

  void AddIdentityMatrix(int n) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        input_data_.push_back(i == j ? 1.0 : 0.0);
      }
    }
  }

  void AddSequentialMatrix(int n) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        input_data_.push_back(static_cast<double>((i * n) + j + 1));
      }
    }
  }

  void AddWeightedMatrix(int n, const std::function<double(int, int)> &weight_func) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        input_data_.push_back(weight_func(i, j));
      }
    }
  }

  bool ValidateMultiplicationResult(OutType &output_data) {
    int n = static_cast<int>(input_data_[0]);

    if (n <= 0) {
      return false;
    }

    size_t expected_size = 1 + (2 * static_cast<size_t>(n) * static_cast<size_t>(n));
    if (input_data_.size() != expected_size) {
      return false;
    }

    Matrix a = ExtractMatrixA(n);
    Matrix b = ExtractMatrixB(n);
    Matrix expected = ComputeExpectedResult(a, b);

    if (output_data.empty() || static_cast<int>(output_data[0]) != n) {
      return false;
    }

    size_t output_expected_size = 1 + (static_cast<size_t>(n) * static_cast<size_t>(n));
    if (output_data.size() != output_expected_size) {
      return false;
    }

    return CompareMatrices(output_data, expected, n);
  }

  [[nodiscard]] Matrix ExtractMatrixA(int n) const {
    Matrix a(n);
    int idx = 1;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        a(i, j) = input_data_[idx++];
      }
    }
    return a;
  }

  [[nodiscard]] Matrix ExtractMatrixB(int n) const {
    Matrix b(n);
    int idx = 1 + (n * n);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        b(i, j) = input_data_[idx++];
      }
    }
    return b;
  }

  static Matrix ComputeExpectedResult(const Matrix &a, const Matrix &b) {
    int n = a.size;
    Matrix expected(n);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        double sum = 0.0;
        for (int k = 0; k < n; ++k) {
          sum += a(i, k) * b(k, j);
        }
        expected(i, j) = sum;
      }
    }
    return expected;
  }

  static bool CompareMatrices(const OutType &output_data, const Matrix &expected, int n) {
    const double eps = 1e-6;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        double result_val = output_data[1 + (i * n) + j];
        double expected_val = expected(i, j);
        if (std::abs(result_val - expected_val) > eps) {
          return false;
        }
      }
    }
    return true;
  }

  InType input_data_;
  int test_number_{0};
};

namespace {

TEST_P(MorozovaSStrassenMultiplicationFuncTests, MatrixMultiplication) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 7> kTestParam = {std::make_tuple(1, "2x2"),         std::make_tuple(2, "4x4"),
                                            std::make_tuple(3, "8x8"),         std::make_tuple(4, "16x16"),
                                            std::make_tuple(5, "32x32"),       std::make_tuple(6, "empty"),
                                            std::make_tuple(7, "invalid_size")};

const auto kTestTasksSEQ = ppc::util::AddFuncTask<MorozovaSStrassenMultiplicationSEQ, InType>(
    kTestParam, PPC_SETTINGS_morozova_s_strassen_multiplication);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksSEQ);
const auto kPerfTestName =
    MorozovaSStrassenMultiplicationFuncTests::PrintFuncTestName<MorozovaSStrassenMultiplicationFuncTests>;

INSTANTIATE_TEST_SUITE_P(StrassenMultiplicationTests, MorozovaSStrassenMultiplicationFuncTests, kGtestValues,
                         kPerfTestName);

}  // namespace

}  // namespace morozova_s_strassen_multiplication
