#include <gtest/gtest.h>

#include "zaharov_g_linear_contrast_stretch/all/include/ops_all.hpp"
#include "zaharov_g_linear_contrast_stretch/common/include/common.hpp"
#include "zaharov_g_linear_contrast_stretch/omp/include/ops_omp.hpp"
#include "zaharov_g_linear_contrast_stretch/seq/include/ops_seq.hpp"
#include "zaharov_g_linear_contrast_stretch/stl/include/ops_stl.hpp"
#include "zaharov_g_linear_contrast_stretch/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace zaharov_g_linear_contrast_stretch {

class ExampleRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 200;
  InType input_data_{};

  void SetUp() override {
    input_data_ = kCount_;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return input_data_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(ExampleRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ZaharovGLinContrStrALL, ZaharovGLinContrStrOMP, ZaharovGLinContrStrSEQ,
                                ZaharovGLinContrStrSTL, ZaharovGLinContrStrTBB>(PPC_SETTINGS_zaharov_g_linear_contrast_stretch);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ExampleRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ExampleRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace zaharov_g_linear_contrast_stretch
