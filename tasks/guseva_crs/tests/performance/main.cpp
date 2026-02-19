#include <gtest/gtest.h>

#include "guseva_crs/common/include/common.hpp"
#include "guseva_crs/common/include/test_reader.hpp"
#include "guseva_crs/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace guseva_crs {

class GusevaMatMulCRSPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;
  OutType output_data_;

  void SetUp() override {
    const auto &filename = ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_guseva_crs), "perf.txt");
    const auto &[a, b, c] = ReadTestFromFile(filename);
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
