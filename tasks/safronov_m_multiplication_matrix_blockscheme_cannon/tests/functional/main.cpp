#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "safronov_m_multiplication_matrix_blockscheme_cannon/common/include/common.hpp"
#include "safronov_m_multiplication_matrix_blockscheme_cannon/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace safronov_m_multiplication_matrix_blocksscheme_cannon {

class SafronovMMultiplicationMatrixBlockSchemeCannonFuncTests
    : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<0>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::make_tuple(std::get<1>(params), std::get<2>(params), std::get<3>(params));
    res_ = std::get<4>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if ((res_.size() * res_[0].size()) != (output_data.size() * output_data[0].size())) {
      return false;
    }
    for (size_t i = 0; i < res_.size(); i++) {
      for (size_t j = 0; j < res_[0].size(); j++) {
        if (std::abs(res_[i][j] - output_data[i][j]) > 1e-10) {
          return false;
        }
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType res_;
};

namespace {

TEST_P(SafronovMMultiplicationMatrixBlockSchemeCannonFuncTests, MultiplicationMatrixBlockSchemeCannon) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 8> kTestParam = {
    std::make_tuple("a", 1, std::vector<std::vector<double>>{{1.0, 2.0}, {3.0, 4.0}},
                    std::vector<std::vector<double>>{{5.0, 6.0}, {7.0, 8.0}},
                    std::vector<std::vector<double>>{{19.0, 22.0}, {43.0, 50.0}}),

    std::make_tuple("b", 2, std::vector<std::vector<double>>(4, std::vector<double>(4, 1.0)),
                    std::vector<std::vector<double>>(4, std::vector<double>(4, 2.0)),
                    std::vector<std::vector<double>>(4, std::vector<double>(4, 8.0))),

    std::make_tuple(
        "c", 2, std::vector<std::vector<double>>{{1, 1, 1, 1}, {2, 2, 2, 2}, {3, 3, 3, 3}, {4, 4, 4, 4}},
        std::vector<std::vector<double>>{{1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}},
        std::vector<std::vector<double>>{{4, 8, 12, 16}, {8, 16, 24, 32}, {12, 24, 36, 48}, {16, 32, 48, 64}}),

    std::make_tuple("d", 3, std::vector<std::vector<double>>{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
                    std::vector<std::vector<double>>{{9, 8, 7}, {6, 5, 4}, {3, 2, 1}},
                    std::vector<std::vector<double>>{{30, 24, 18}, {84, 69, 54}, {138, 114, 90}}),

    std::make_tuple("e", 2, std::vector<std::vector<double>>(4, std::vector<double>(4, 0.1)),
                    std::vector<std::vector<double>>(4, std::vector<double>(4, 0.2)),
                    std::vector<std::vector<double>>(4, std::vector<double>(4, 0.08))),

    std::make_tuple("f", 3, std::vector<std::vector<double>>(6, std::vector<double>(6, 1.5)),
                    std::vector<std::vector<double>>(6, std::vector<double>(6, 0.5)),
                    std::vector<std::vector<double>>(6, std::vector<double>(6, 4.5))),

    std::make_tuple("g", 4, std::vector<std::vector<double>>(8, std::vector<double>(8, 1.25)),
                    std::vector<std::vector<double>>(8, std::vector<double>(8, 0.8)),
                    std::vector<std::vector<double>>(8, std::vector<double>(8, 8.0))),

    std::make_tuple("h", 3, std::vector<std::vector<double>>(9, std::vector<double>(9, 1.1)),
                    std::vector<std::vector<double>>(9, std::vector<double>(9, 2.0)),
                    std::vector<std::vector<double>>(9, std::vector<double>(9, 19.8)))};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<SafronovMMultiplicationMatrixBlockSchemeCannon, InType>(
        kTestParam, PPC_SETTINGS_safronov_m_multiplication_matrix_blockscheme_cannon));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = SafronovMMultiplicationMatrixBlockSchemeCannonFuncTests::PrintFuncTestName<
    SafronovMMultiplicationMatrixBlockSchemeCannonFuncTests>;

INSTANTIATE_TEST_SUITE_P(MultiplicationMatrixBlockSchemeCannonTests,
                         SafronovMMultiplicationMatrixBlockSchemeCannonFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace safronov_m_multiplication_matrix_blocksscheme_cannon
