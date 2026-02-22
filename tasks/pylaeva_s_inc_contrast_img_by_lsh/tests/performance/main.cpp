#include <gtest/gtest.h>

#include "pylaeva_s_inc_contrast_img_by_lsh/all/include/ops_all.hpp"
#include "pylaeva_s_inc_contrast_img_by_lsh/common/include/common.hpp"
#include "pylaeva_s_inc_contrast_img_by_lsh/omp/include/ops_omp.hpp"
#include "pylaeva_s_inc_contrast_img_by_lsh/seq/include/ops_seq.hpp"
#include "pylaeva_s_inc_contrast_img_by_lsh/stl/include/ops_stl.hpp"
#include "pylaeva_s_inc_contrast_img_by_lsh/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace pylaeva_s_inc_contrast_img_by_lsh {

class PylaevaSRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
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

TEST_P(PylaevaSRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, PylaevaSIncContrastImgByLshALL, PylaevaSIncContrastImgByLshOMP, PylaevaSIncContrastImgByLshSEQ,
                                PylaevaSIncContrastImgByLshSTL, PylaevaSIncContrastImgByLshTBB>(PPC_SETTINGS_pylaeva_s_inc_contrast_img_by_lsh);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = PylaevaSRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, PylaevaSRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace pylaeva_s_inc_contrast_img_by_lsh
