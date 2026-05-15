#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "baranov_a_mult_matrix_fox_algorithm/common/include/common.hpp"
#include "baranov_a_mult_matrix_fox_algorithm/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace baranov_a_mult_matrix_fox_algorithm_seq {

class BaranovAFuncTest : public ppc::util::BaseRunFuncTests<baranov_a_mult_matrix_fox_algorithm::InType,
                                                            baranov_a_mult_matrix_fox_algorithm::OutType,
                                                            baranov_a_mult_matrix_fox_algorithm::TestType> {
 public:
  static std::string PrintTestParam(const baranov_a_mult_matrix_fox_algorithm::TestType &test_param) {
    size_t n = std::get<0>(test_param);
    std::string type = std::get<1>(test_param);
    return "n_" + std::to_string(n) + "_" + type;
  }

 protected:
  void SetUp() override {
    auto params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    size_t n = std::get<0>(params);
    std::string test_type = std::get<1>(params);
    size_t size = n * n;

    std::vector<double> a(size);
    std::vector<double> b(size);

    if (test_type.find("identity") != std::string::npos) {
      GenerateIdentityMatrix(a, b, n);
    } else if (test_type.find("random") != std::string::npos) {
      GenerateRandomMatrices(a, b, n, std::stoi(test_type.substr(test_type.find("_seed") + 5)));
    } else if (test_type.find("extreme") != std::string::npos) {
      GenerateExtremeValuesMatrices(a, b, n);
    } else if (test_type.find("sparse") != std::string::npos) {
      GenerateSparseMatrices(a, b, n);
    } else if (test_type.find("constant") != std::string::npos) {
      GenerateConstantMatrices(a, b, n);
    } else {
      GenerateLinearMatrices(a, b, n);
    }

    input_data_ = std::make_tuple(n, a, b);
    std::vector<double> expected(size, 0.0);
    ReferenceMultiply(a, b, expected, n);
    expected_output_ = expected;
  }

  bool CheckTestOutputData(baranov_a_mult_matrix_fox_algorithm::OutType &output_data) final {
    if (expected_output_.size() != output_data.size()) {
      ADD_FAILURE() << "Size mismatch: expected " << expected_output_.size() << ", got " << output_data.size();
      return false;
    }

    const double epsilon = 1e-8;
    size_t mismatches = 0;
    double max_diff = 0.0;

    for (size_t i = 0; i < expected_output_.size(); ++i) {
      double diff = std::abs(expected_output_[i] - output_data[i]);
      max_diff = std::max(max_diff, diff);
      if (diff > epsilon) {
        mismatches++;
        if (mismatches <= 5) {
          ADD_FAILURE() << "Mismatch at index " << i << ": expected " << expected_output_[i] << ", got "
                        << output_data[i];
        }
      }
    }

    if (mismatches > 0) {
      ADD_FAILURE() << "Total mismatches: " << mismatches << ", max difference: " << max_diff;
      return false;
    }
    return true;
  }

  baranov_a_mult_matrix_fox_algorithm::InType GetTestInputData() final {
    return input_data_;
  }

 private:
  static void ReferenceMultiply(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c,
                                size_t n) {
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        double sum = 0.0;
        for (size_t k = 0; k < n; ++k) {
          sum += a[(i * n) + k] * b[(k * n) + j];
        }
        c[(i * n) + j] = sum;
      }
    }
  }

  static void GenerateLinearMatrices(std::vector<double> &a, std::vector<double> &b, size_t n) {
    size_t size = n * n;
    for (size_t i = 0; i < size; ++i) {
      a[i] = static_cast<double>(i + 1);
      b[i] = static_cast<double>(size - i);
    }
  }

  static void GenerateIdentityMatrix(std::vector<double> &a, std::vector<double> &b, size_t n) {
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        a[(i * n) + j] = (i == j) ? 1.0 : 0.0;
      }
    }
    for (size_t i = 0; i < n * n; ++i) {
      b[i] = static_cast<double>(i + 1) * 0.5;
    }
  }

  static void GenerateRandomMatrices(std::vector<double> &a, std::vector<double> &b, size_t n, int seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-100.0, 100.0);

    for (size_t i = 0; i < n * n; ++i) {
      a[i] = dist(gen);
      b[i] = dist(gen);
    }
  }

  static void GenerateExtremeValuesMatrices(std::vector<double> &a, std::vector<double> &b, size_t n) {
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        if ((i + j) % 2 == 0) {
          a[(i * n) + j] = 1e6 * static_cast<double>(i + 1);
          b[(i * n) + j] = 1e-6 * static_cast<double>(j + 1);
        } else {
          a[(i * n) + j] = 1e-6 * static_cast<double>(i + 1);
          b[(i * n) + j] = 1e6 * static_cast<double>(j + 1);
        }
      }
    }
  }

  static void GenerateSparseMatrices(std::vector<double> &a, std::vector<double> &b, size_t n) {
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        a[(i * n) + j] = 0.0;
        b[(i * n) + j] = 0.0;
      }
    }

    for (size_t i = 0; i < n; ++i) {
      a[(i * n) + i] = static_cast<double>(i + 1);
      b[(i * n) + i] = static_cast<double>(n - i);
    }

    if (n > 2) {
      a[(0 * n) + 2] = 5.0;
      a[(1 * n) + 3] = 7.0;
      b[(2 * n) + 1] = 3.0;
      b[(3 * n) + 0] = 4.0;
    }
  }

  static void GenerateConstantMatrices(std::vector<double> &a, std::vector<double> &b, size_t n) {
    double const_a = 2.5;
    double const_b = 1.5;

    for (size_t i = 0; i < n * n; ++i) {
      a[i] = const_a;
      b[i] = const_b;
    }
  }

  baranov_a_mult_matrix_fox_algorithm::InType input_data_;
  baranov_a_mult_matrix_fox_algorithm::OutType expected_output_;
};

namespace {

TEST_P(BaranovAFuncTest, MatrixMultiplicationTest) {
  ExecuteTest(GetParam());
}

const std::array<baranov_a_mult_matrix_fox_algorithm::TestType, 20> kTestParams = {
    std::make_tuple(1, "size1_simple"),   std::make_tuple(2, "size2_simple"),   std::make_tuple(3, "size3_simple"),

    std::make_tuple(2, "identity_2"),     std::make_tuple(4, "identity_4"),     std::make_tuple(8, "identity_8"),

    std::make_tuple(3, "random_seed123"), std::make_tuple(5, "random_seed456"), std::make_tuple(7, "random_seed789"),

    std::make_tuple(4, "extreme_4"),      std::make_tuple(6, "extreme_6"),

    std::make_tuple(4, "sparse_4"),       std::make_tuple(8, "sparse_8"),

    std::make_tuple(3, "constant_3"),     std::make_tuple(5, "constant_5"),     std::make_tuple(7, "constant_7"),

    std::make_tuple(16, "size16_block"),  std::make_tuple(32, "size32_block"),  std::make_tuple(64, "size64_block"),
    std::make_tuple(128, "size128_block")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<BaranovAMultMatrixFoxAlgorithmSEQ, baranov_a_mult_matrix_fox_algorithm::InType>(
        kTestParams, PPC_SETTINGS_baranov_a_mult_matrix_fox_algorithm);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = BaranovAFuncTest::PrintFuncTestName<BaranovAFuncTest>;

INSTANTIATE_TEST_SUITE_P(FoxAlgorithmTests, BaranovAFuncTest, kGtestValues, kTestName);

}  // namespace

}  // namespace baranov_a_mult_matrix_fox_algorithm_seq
