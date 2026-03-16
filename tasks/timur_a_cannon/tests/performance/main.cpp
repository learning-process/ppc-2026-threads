#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "timur_a_cannon/common/include/common.hpp"
#include "timur_a_cannon/omp/include/ops_omp.hpp"
#include "timur_a_cannon/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace timur_a_cannon {

class TimurACannonPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  InType input_data_;

  void SetUp() override {
    const int kCount = 512;
    int size_block = 32;
    std::vector<std::vector<double>> matrix_a(kCount, std::vector<double>(kCount, 2.0));
    std::vector<std::vector<double>> matrix_b(kCount, std::vector<double>(kCount, 3.0));

    input_data_ = std::make_tuple(size_block, matrix_a, matrix_b);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return !output_data.empty() && !output_data[0].empty();
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(TimurACannonPerfTests, MultiplicationMatrixBlockSchemeCannonPerf) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, TimurACannonMatrixMultiplication, TimurACannonMatrixMultiplicationOMP>(
        PPC_SETTINGS_timur_a_cannon);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = TimurACannonPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, TimurACannonPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace timur_a_cannon
