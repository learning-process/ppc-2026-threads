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

#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"
#include "zagryadskov_m_complex_spmm_ccs/common/include/common.hpp"
#include "zagryadskov_m_complex_spmm_ccs/seq/include/ops_seq.hpp"

namespace zagryadskov_m_complex_spmm_ccs {

class ZagryadskovMRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
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
      a.values = {1.0, 2.0, 3.0};

      b.m = 3;
      b.n = 2;
      b.col_ptr = {0, 2, 3};
      b.row_ind = {0, 2, 1};
      b.values = {4.0, 5.0, 6.0};

      c.m = 2;
      c.n = 2;
      c.col_ptr = {0, 1, 2};
      c.row_ind = {0, 1};
      c.values = {19.0, 12.0};
    }

    if (params == 1) {
      a.m = 2;
      a.n = 3;
      a.col_ptr = {0, 1, 2, 4};
      a.row_ind = {0, 1, 0, 1};
      a.values = {1.0, 3.0, 2.0, 4.0};

      b.m = 3;
      b.n = 2;
      b.col_ptr = {0, 2, 4};
      b.row_ind = {0, 1, 1, 2};
      b.values = {5.0, 6.0, 7.0, 8.0};

      c.m = 2;
      c.n = 2;
      c.col_ptr = {0, 2, 4};
      c.row_ind = {0, 1, 0, 1};
      c.values = {5.0, 18.0, 16.0, 53.0};
    }

    if (params == 2) {
      a.m = 3;
      a.n = 3;
      a.col_ptr = {0, 1, 2, 3};
      a.row_ind = {0, 1, 2};
      a.values = {1.0, 2.0, 3.0};

      b.m = 3;
      b.n = 3;
      b.col_ptr = {0, 1, 2, 3};
      b.row_ind = {2, 0, 2};
      b.values = {5.0, 4.0, 6.0};

      c.m = 3;
      c.n = 3;
      c.col_ptr = {0, 1, 2, 3};
      c.row_ind = {2, 0, 2};
      c.values = {15.0, 4.0, 18.0};
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    bool result = true;
    double eps = 1e-14;
    bool f1 = test_result_.m != output_data.m;
    bool f2 = test_result_.n != output_data.n;
    bool f3 = test_result_.col_ptr.size() != output_data.col_ptr.size();
    bool f4 = test_result_.row_ind.size() != output_data.row_ind.size();
    bool f5 = test_result_.values.size() != output_data.values.size();

    if (f1 || f2 || f3 || f4 || f5) {
      result = false;
    }
    for (size_t i = 0; i < test_result_.col_ptr.size(); ++i) {
      if (test_result_.col_ptr[i] != output_data.col_ptr[i]) {
        result = false;
      }
    }

    for (int j = 0; j < test_result_.n; ++j) {
      std::vector<std::pair<int, std::complex<double>>> test;
      std::vector<std::pair<int, std::complex<double>>> output;
      for (int k = test_result_.col_ptr[j]; k < test_result_.col_ptr[j + 1]; ++k) {
        test.emplace_back(test_result_.row_ind[k], test_result_.values[k]);
        output.emplace_back(output_data.row_ind[k], output_data.values[k]);
      }
      auto cmp = [](const auto &x, const auto &y) { return x.first < y.first; };
      std::ranges::sort(test, cmp);
      std::ranges::sort(output, cmp);

      for (size_t i = 0; i < test.size(); ++i) {
        bool f6 = test[i].first != output[i].first;
        bool f7 = std::abs(test[i].second - output[i].second) > eps;
        if (f6 || f7) {
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

TEST_P(ZagryadskovMRunFuncTestsThreads, FuncCCSTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {0, 1, 2};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<ZagryadskovMComplexSpMMCCSSEQ, InType>(
    kTestParam, PPC_SETTINGS_zagryadskov_m_complex_spmm_ccs));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = ZagryadskovMRunFuncTestsThreads::PrintFuncTestName<ZagryadskovMRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(RunFuncCCSTest, ZagryadskovMRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace zagryadskov_m_complex_spmm_ccs
