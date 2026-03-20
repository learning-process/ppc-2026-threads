#include <gtest/gtest.h>

#include <cmath>

#include "ovsyannikov_n_simpson_method/common/include/common.hpp"
#include "ovsyannikov_n_simpson_method/omp/include/ops_omp.hpp"
#include "ovsyannikov_n_simpson_method/seq/include/ops_seq.hpp"
#include "ovsyannikov_n_simpson_method/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace ovsyannikov_n_simpson_method {

class OvsyannikovNRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  void SetUp() override {
    input_data_ = InType{0.0, 1.0, 0.0, 1.0, 2000, 2000};
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::abs(output_data - 1.0) < 1e-4;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = {};
};

TEST_P(OvsyannikovNRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kPerfTasksSEQ =
    ppc::util::MakeAllPerfTasks<InType, OvsyannikovNSimpsonMethodSEQ>(PPC_SETTINGS_ovsyannikov_n_simpson_method);

const auto kPerfTasksOMP =
    ppc::util::MakeAllPerfTasks<InType, OvsyannikovNSimpsonMethodOMP>(PPC_SETTINGS_ovsyannikov_n_simpson_method);

const auto kPerfTasksTBB =
    ppc::util::MakeAllPerfTasks<InType, OvsyannikovNSimpsonMethodTBB>(PPC_SETTINGS_ovsyannikov_n_simpson_method);

const auto kPerfTestName = OvsyannikovNRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(SimpsonPerf_SEQ, OvsyannikovNRunPerfTestThreads, ppc::util::TupleToGTestValues(kPerfTasksSEQ),
                         kPerfTestName);

INSTANTIATE_TEST_SUITE_P(SimpsonPerf_OMP, OvsyannikovNRunPerfTestThreads, ppc::util::TupleToGTestValues(kPerfTasksOMP),
                         kPerfTestName);

INSTANTIATE_TEST_SUITE_P(SimpsonPerf_TBB, OvsyannikovNRunPerfTestThreads, ppc::util::TupleToGTestValues(kPerfTasksTBB),
                         kPerfTestName);

}  // namespace
}  // namespace ovsyannikov_n_simpson_method
