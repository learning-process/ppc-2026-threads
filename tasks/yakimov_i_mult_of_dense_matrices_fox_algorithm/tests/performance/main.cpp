#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>

#include "util/include/perf_test_util.hpp"
#include "yakimov_i_mult_of_dense_matrices_fox_algorithm_seq/common/include/common.hpp"
#include "yakimov_i_mult_of_dense_matrices_fox_algorithm_seq/seq/include/ops_seq.hpp"

namespace yakimov_i_mult_of_dense_matrices_fox_algorithm_seq {

class YakimovIMultDenseFoxPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  bool CheckTestOutputData(OutType &output_data) final {
    return std::isfinite(output_data);
  }

  InType GetTestInputData() final {
    static size_t test_index = 0;
    static constexpr std::array<InType, 4> kTestSizes = {10, 16, 20, 32};
    InType result = kTestSizes.at(test_index % kTestSizes.size());
    test_index++;
    return result;
  }
};

TEST_P(YakimovIMultDenseFoxPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, YakimovIMultOfDenseMatricesFoxAlgorithmSEQ>(
    PPC_SETTINGS_yakimov_i_mult_of_dense_matrices_fox_algorithm_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = YakimovIMultDenseFoxPerfTests::CustomPerfTestName;

namespace {

INSTANTIATE_TEST_SUITE_P(RunModeTests, YakimovIMultDenseFoxPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace yakimov_i_mult_of_dense_matrices_fox_algorithm_seq
