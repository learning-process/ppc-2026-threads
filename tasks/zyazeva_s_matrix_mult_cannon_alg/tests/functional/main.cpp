#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"
#include "zyazeva_s_matrix_mult_cannon_alg/common/include/common.hpp"
#include "zyazeva_s_matrix_mult_cannon_alg/seq/include/ops_seq.hpp"

namespace zyazeva_s_matrix_mult_cannon_alg {

class ZyazevaARunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param) + "_" + std::to_string(std::get<0>(test_param)) + "by" +
           std::to_string(std::get<0>(test_param)) + "_" + std::to_string(std::get<1>(test_param));
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    size_t sz = std::get<0>(params);
    size_t size = sz * sz;

    int up_to = std::get<1>(params);

    std::vector<double> m1(sz * sz);
    std::vector<double> m2(sz * sz);

    for (size_t i = 0; i < size; i++) {
      m1[i] = static_cast<double>(up_to - static_cast<int>(i));
      m2[i] = static_cast<double>(i) * static_cast<double>(up_to) / static_cast<double>(size - 1);
    }

    input_data_ = std::make_tuple(sz, m1, m2);

    std::vector<double> res(size, 0.0);
    MatMul(m1, m2, res, static_cast<int>(sz));
    expected_output_ = res;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (expected_output_.size() != output_data.size()) {
      return false;
    }
    const double epsilon = 1e-7;
    for (size_t i = 0; i < expected_output_.size(); ++i) {
      if (std::abs(expected_output_[i] - output_data[i]) > epsilon) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

  static void MatMul(const std::vector<double> &m1, const std::vector<double> &m2, std::vector<double> &res, int n) {
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < n; ++k) {
        double tmp = m1[(i * n) + k];
        for (int j = 0; j < n; ++j) {
          res[(i * n) + j] += tmp * m2[(k * n) + j];
        }
      }
    }
  }

 private:
  InType input_data_;
  OutType expected_output_;
};

namespace {

TEST_P(ZyazevaARunFuncTestsThreads, MatMulCannonAlg) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 7> kTestParams = {
    std::make_tuple(1, 70, "small_matrix"),        std::make_tuple(2, 100, "small2_matrix"),
    std::make_tuple(3, 120, "medium1_matrix"),     std::make_tuple(4, 150, "medium2_matrix"),
    std::make_tuple(7, 300, "medium3_matrix"),     std::make_tuple(10, 700, "large_matrix"),
    std::make_tuple(15, 1500, "very_large_matrix")};

const auto kTestTasksList = ppc::util::AddFuncTask<ZyazevaSMatrixMultCannonAlgSEQ, InType>(
    kTestParams, PPC_SETTINGS_zyazeva_s_matrix_mult_cannon_alg);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = ZyazevaARunFuncTestsThreads::PrintFuncTestName<ZyazevaARunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(MatMulCannonAlg, ZyazevaARunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace zyazeva_s_matrix_mult_cannon_alg
