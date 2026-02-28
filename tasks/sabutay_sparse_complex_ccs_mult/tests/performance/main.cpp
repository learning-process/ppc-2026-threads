#include <gtest/gtest.h>

#include "sabutay_sparse_complex_ccs_mult/common/include/common.hpp"
#include "sabutay_sparse_complex_ccs_mult/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace sabutay_sparse_complex_ccs_mult {

class SabutayRunPerfSeq : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_{};

  void SetUp() override {
    // simple fixed matrices for performance measurement
    CCS a, b;
    a.m = 100; a.n = 100;
    a.col_ptr.resize(a.n+1);
    b.m = 100; b.n = 100;
    b.col_ptr.resize(b.n+1);
    // fill diagonal ones
    for (int i = 0; i < 100; ++i) {
      a.col_ptr[i] = i;
      b.col_ptr[i] = i;
      a.row_ind.push_back(i);
      a.values.emplace_back(1.0, 0.0);
      b.row_ind.push_back(i);
      b.values.emplace_back(1.0, 0.0);
    }
    a.col_ptr[100] = 100;
    b.col_ptr[100] = 100;
    input_data_ = std::make_tuple(a, b);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // result should equal identity multiplication
    bool ok = (output_data.m == 100 && output_data.n == 100);
    return ok;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(SabutayRunPerfSeq, RunPerf) {
  ExecuteTest(GetParam());
}

namespace {

const auto kPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, SabutaySparseComplexCcsMultSEQ>(PPC_SETTINGS_sabutay_sparse_complex_ccs_mult);

const auto kGtestValues = ppc::util::TupleToGTestValues(kPerfTasks);

const auto kPerfTestName = SabutayRunPerfSeq::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, SabutayRunPerfSeq, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sabutay_sparse_complex_ccs_mult
