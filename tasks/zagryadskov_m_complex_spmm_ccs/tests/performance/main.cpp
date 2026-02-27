#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <random>
#include <utility>
#include <vector>

#include "util/include/perf_test_util.hpp"
#include "zagryadskov_m_complex_spmm_ccs/common/include/common.hpp"
#include "zagryadskov_m_complex_spmm_ccs/seq/include/ops_seq.hpp"

namespace zagryadskov_m_complex_spmm_ccs {

class ZagryadskovMRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;
  OutType test_result_;

  void SetUp() override {
    int dim = 20000;
    int seed = 0;
    CCS &a = std::get<0>(input_data_);
    CCS &b = std::get<1>(input_data_);
    CCS &c = test_result_;

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> val_gen(1.0, 2.0);
    std::uniform_int_distribution<int> size_gen(0, 100);
    std::vector<int> indices(dim);
    for (size_t i = 0; i < indices.size(); ++i) {
      indices[i] = static_cast<int>(i);
    }

    a.m = dim;
    a.n = dim;
    a.col_ptr.assign(a.n + 1, 0);
    for (int j = 0; j < a.n; ++j) {
      int size = size_gen(rng);
      a.col_ptr[j + 1] = a.col_ptr[j] + size;
      std::shuffle(indices.begin(), indices.end(), rng);
      for (int k = 0; k < size; ++k) {
        int ind = indices[k];
        double av = val_gen(rng);
        double bv = val_gen(rng);
        std::complex<double> z(av, bv);
        a.row_ind.push_back(ind);
        a.values.push_back(z);
      }
    }

    b.m = dim;
    b.n = dim;
    b.col_ptr.assign(b.n + 1, 0);
    for (int j = 0; j < b.n; ++j) {
      int size = size_gen(rng);
      b.col_ptr[j + 1] = b.col_ptr[j] + size;
      std::shuffle(indices.begin(), indices.end(), rng);
      for (int k = 0; k < size; ++k) {
        int ind = indices[k];
        double av = val_gen(rng);
        double bv = val_gen(rng);
        std::complex<double> z(av, bv);
        b.row_ind.push_back(ind);
        b.values.push_back(z);
      }
    }

    ZagryadskovMComplexSpMMCCSSEQ::SpMM(a, b, c);
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
};

TEST_P(ZagryadskovMRunPerfTestThreads, PerfCCSTest) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ZagryadskovMComplexSpMMCCSSEQ>(PPC_SETTINGS_zagryadskov_m_complex_spmm_ccs);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ZagryadskovMRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunPerfCCSTest, ZagryadskovMRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace zagryadskov_m_complex_spmm_ccs
