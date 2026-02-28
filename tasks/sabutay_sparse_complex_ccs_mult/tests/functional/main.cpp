#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "sabutay_sparse_complex_ccs_mult/common/include/common.hpp"
#include "sabutay_sparse_complex_ccs_mult/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace sabutay_sparse_complex_ccs_mult {

class SabutayARunFuncTestsSeq : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) { return std::to_string(test_param); }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    CCS &a = std::get<0>(input_data_);
    CCS &b = std::get<1>(input_data_);
    CCS &c = test_result_;

    if (params == 0) {
      a.m = 2; a.n = 3;
      a.col_ptr = {0,1,2,3};
      a.row_ind = {0,1,0};
      a.values = {{1.0,0.0},{2.0,0.0},{3.0,0.0}};

      b.m = 3; b.n = 2;
      b.col_ptr = {0,2,3};
      b.row_ind = {0,2,1};
      b.values = {{4.0,0.0},{5.0,0.0},{6.0,0.0}};

      c.m = 2; c.n = 2;
      c.col_ptr = {0,1,2};
      c.row_ind = {0,1};
      c.values = {{19.0,0.0},{12.0,0.0}};
    }
    if (params == 1) {
      a.m = 2; a.n = 2;
      a.col_ptr = {0,2,3};
      a.row_ind = {0,1,1};
      a.values = {{1.0,1.0},{-1.0,2.0},{3.0,0.0}};

      b.m = 2; b.n = 2;
      b.col_ptr = {0,1,3};
      b.row_ind = {1,0,1};
      b.values = {{2.0,0.5},{-1.0,1.0},{0.0,-1.0}};

      c.m = 2; c.n = 2;
      c.col_ptr = {0,1,3};
      c.row_ind = {1,0,1};
      c.values = {{-2.0,2.5},{-5.0,1.0},{-3.0,5.0}};
    }
    if (params == 2) {
      a.m = 1; a.n = 1;
      a.col_ptr = {0,1};
      a.row_ind = {0};
      a.values = {{3.0,-1.0}};

      b.m = 1; b.n = 1;
      b.col_ptr = {0,1};
      b.row_ind = {0};
      b.values = {{0.0,2.0}};

      c.m = 1; c.n = 1;
      c.col_ptr = {0,1};
      c.row_ind = {0};
      c.values = {{2.0,6.0}};
    }
  }

  bool CheckTestOutputData(OutType& output_data) final {
    bool result = true;
    double eps = 1e-14;
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
      for (int k = test_result_.col_ptr[j]; k < test_result_.col_ptr[j+1]; ++k) {
        test.emplace_back(test_result_.row_ind[k], test_result_.values[k]);
        output.emplace_back(output_data.row_ind[k], output_data.values[k]);
      }
      auto cmp = [](auto const &x, auto const &y){ return x.first < y.first; };
      std::ranges::sort(test, cmp);
      std::ranges::sort(output, cmp);
      for (size_t i = 0; i < test.size(); ++i) {
        if (test[i].first != output[i].first ||
            std::abs(test[i].second - output[i].second) > eps) {
          result = false;
        }
      }
    }

    return result;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType test_result_;
};

namespace {

TEST_P(SabutayARunFuncTestsSeq, FuncCCSTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {0,1,2};

const auto kTestTasksList = ppc::util::AddFuncTask<SabutaySparseComplexCcsMultSEQ, InType>(kTestParam, PPC_SETTINGS_sabutay_sparse_complex_ccs_mult);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = SabutayARunFuncTestsSeq::PrintFuncTestName<SabutayARunFuncTestsSeq>;

INSTANTIATE_TEST_SUITE_P(RunFuncCCSTest, SabutayARunFuncTestsSeq, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sabutay_sparse_complex_ccs_mult