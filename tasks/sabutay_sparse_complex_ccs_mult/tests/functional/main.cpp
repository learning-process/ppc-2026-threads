#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <string>

#include "sabutay_sparse_complex_ccs_mult/common/include/common.hpp"
#include "sabutay_sparse_complex_ccs_mult/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace sabutay_sparse_complex_ccs_mult {

// Constants for test values 
constexpr double kZero = 0.0;
constexpr double kValue1 = 1.0;
constexpr double kValue2 = 2.0;
constexpr double kValue3 = 3.0;
constexpr double kValue4 = 4.0;
constexpr double kValue5 = 5.0;
constexpr double kValue6 = 6.0;
constexpr double kValue12 = 12.0;
constexpr double kValue19 = 19.0;

class SabutayARunFuncTestsSeq : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static auto PrintTestParam(const TestType &test_param) -> std::string {
    return std::to_string(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    CCS &a = std::get<0>(input_data_);
    CCS &b = std::get<1>(input_data_);
    CCS &c = test_result_;

    if (params == 0) {
      a.m = 2;
      a.n = 3;
      a.col_ptr = {0, 1, 2, 3};
      a.row_ind = {0, 1, 0};
      a.values = {{kValue1, kZero}, {kValue2, kZero}, {kValue3, kZero}};

      b.m = 3;
      b.n = 2;
      b.col_ptr = {0, 2, 3};
      b.row_ind = {0, 2, 1};
      b.values = {{kValue4, kZero}, {kValue5, kZero}, {kValue6, kZero}};

      c.m = 2;
      c.n = 2;
      c.col_ptr = {0, 1, 2};
      c.row_ind = {0, 1};
      c.values = {{kValue19, kZero}, {kValue12, kZero}};
    }
    if (params == 1) {
      a.m = 2;
      a.n = 2;
      a.col_ptr = {0, 1, 2};
      a.row_ind = {0, 1};
      a.values = {{kValue1, kZero}, {kValue2, kZero}};

      b.m = 2;
      b.n = 2;
      b.col_ptr = {0, 1, 2};
      b.row_ind = {1, 0};
      b.values = {{kValue3, kZero}, {kValue4, kZero}};

      c.m = 2;
      c.n = 2;
      c.col_ptr = {0, 1, 2};
      c.row_ind = {1, 0};
      c.values = {{kValue6, kZero}, {kValue4, kZero}};
    }
    if (params == 2) {
      a.m = 1;
      a.n = 1;
      a.col_ptr = {0, 1};
      a.row_ind = {0};
      a.values = {{kValue2, kValue1}};

      b.m = 1;
      b.n = 1;
      b.col_ptr = {0, 1};
      b.row_ind = {0};
      b.values = {{kValue1, kValue1}};

      c.m = 1;
      c.n = 1;
      c.col_ptr = {0, 1};
      c.row_ind = {0};
      c.values = {{kValue1, kValue3}};
    }
  }

  bool CheckTestOutputData(OutType &output_data) {
    bool result = true;
    constexpr double kEps = 1e-14;
    if (test_result_.m != output_data.m || test_result_.n != output_data.n ||
        test_result_.col_ptr.size() != output_data.col_ptr.size() ||
        test_result_.row_ind.size() != output_data.row_ind.size() ||
        test_result_.values.size() != output_data.values.size()) {
      return false;
    }

    for (size_t i = 0; i < test_result_.col_ptr.size(); ++i) {
      if (test_result_.col_ptr[i] != output_data.col_ptr[i]) {
        return false;
      }
    }

    for (int j = 0; j < test_result_.n; ++j) {
      std::vector<std::pair<int, std::complex<double>>> test;
      std::vector<std::pair<int, std::complex<double>>> output;
      for (int k = 0; k < test_result_.col_ptr[j + 1] - test_result_.col_ptr[j]; ++k) {
        test.emplace_back(test_result_.row_ind[test_result_.col_ptr[j] + k],
                          test_result_.values[test_result_.col_ptr[j] + k]);
        output.emplace_back(output_data.row_ind[test_result_.col_ptr[j] + k],
                            output_data.values[test_result_.col_ptr[j] + k]);
      }
      auto cmp = [](const auto &x, const auto &y) { return x.first < y.first; };
      std::sort(test.begin(), test.end(), cmp);
      std::sort(output.begin(), output.end(), cmp);
      for (size_t i = 0; i < test.size(); ++i) {
        if (test[i].first != output[i].first || std::abs(test[i].second - output[i].second) > kEps) {
          result = false;
        }
      }
    }

    return result;
  }

  InType GetTestInputData() {
    return input_data_;
  }

 private:
  InType input_data_{};
  OutType test_result_{};
};

namespace {

TEST_P(SabutayARunFuncTestsSeq, FuncCCSTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {0, 1, 2};

const auto kTestTasksList = ppc::util::AddFuncTask<SabutaySparseComplexCcsMultSEQ, InType>(
    kTestParam, PPC_SETTINGS_sabutay_sparse_complex_ccs_mult);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = SabutayARunFuncTestsSeq::PrintFuncTestName<SabutayARunFuncTestsSeq>;

INSTANTIATE_TEST_SUITE_P(RunFuncCCSTest, SabutayARunFuncTestsSeq, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sabutay_sparse_complex_ccs_mult
