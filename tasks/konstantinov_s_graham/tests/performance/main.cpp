#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

// #include "konstantinov_s_graham/all/include/ops_all.hpp"
#include "konstantinov_s_graham/common/include/common.hpp"
// #include "konstantinov_s_graham/omp/include/ops_omp.hpp"
#include "konstantinov_s_graham/seq/include/ops_seq.hpp"
// #include "konstantinov_s_graham/stl/include/ops_stl.hpp"
// #include "konstantinov_s_graham/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace konstantinov_s_graham {

class KonstantinovSRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType test_input_;
  OutType test_expected_output_;

  void SetUp() override {
    const size_t side = 1500;
    const size_t total = side * side;

    test_input_.first.reserve(total);
    test_input_.second.reserve(total);

    for (size_t i = 0; i < side; ++i) {
      for (size_t j = 0; j < side; ++j) {
        test_input_.first.push_back(static_cast<double>(i));
        test_input_.second.push_back(static_cast<double>(j));
      }
    }

    test_expected_output_ = {{0.0, 0.0},
                             {static_cast<double>(side - 1), 0.0},
                             {static_cast<double>(side - 1), static_cast<double>(side - 1)},
                             {0.0, static_cast<double>(side - 1)}};
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.size() == 4;
  }

  InType GetTestInputData() final {
    return test_input_;
  }
};

TEST_P(KonstantinovSRunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KonstantinovAGrahamSEQ>(PPC_SETTINGS_konstantinov_s_graham);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KonstantinovSRunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KonstantinovSRunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace konstantinov_s_graham
