#include <gtest/gtest.h>

#include <cstddef>
#include <tuple>
#include <vector>

#include "guseva_crs/common/include/common.hpp"
#include "guseva_crs/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace guseva_crs {

class GusevaMatMulCRSPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;
  OutType output_data_;
  std::size_t inp_size_ = 10000;

  static CRS CreateDiagMatrix(double value, std::size_t size) {
    CRS a;
    a.ncols = size;
    a.nrows = size;
    a.nz = size;
    a.values = std::vector<double>(a.nz, value);
    a.cols = std::vector<std::size_t>(a.nz, 0);
    for (std::size_t i = 0; i < a.nz; i++) {
      a.cols[i] = i;
    }
    a.row_ptrs = std::vector<std::size_t>(a.nrows + 1, 0);
    for (std::size_t i = 0; i < a.nrows + 1; i++) {
      a.row_ptrs[i] = i;
    }
    return a;
  }

  void SetUp() override {
    CRS a = CreateDiagMatrix(2.5, inp_size_);
    CRS b = CreateDiagMatrix(2.5, inp_size_);
    CRS c = CreateDiagMatrix(6.25, inp_size_);
    input_data_ = std::make_tuple(a, b);
    output_data_ = c;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return Equal(output_data, output_data_);
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(GusevaMatMulCRSPerfTest, G) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, GusevaCRSMatMulSeq>(PPC_SETTINGS_guseva_crs);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = GusevaMatMulCRSPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, GusevaMatMulCRSPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace guseva_crs
