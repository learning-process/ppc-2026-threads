#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <random>

#include "lopatin_a_sobel_operator/common/include/common.hpp"
#include "lopatin_a_sobel_operator/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace lopatin_a_sobel_operator {

class LopatinARunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;
  OutType output_chekup_data_;

  void SetUp() override {
    std::size_t height = 3840;
    std::size_t width = 2160;

    input_data_.height = height;
    input_data_.width = width;
    input_data_.threshold = 100;
    input_data_.pixels.resize(height * width);

    int seed = 20260223;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<int> dis(0, 255);

    for (unsigned char &pixel : input_data_.pixels) {
      pixel = static_cast<std::uint8_t>(dis(gen));
    }

    output_chekup_data_.resize(height * width);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.size() == output_chekup_data_.size();
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(LopatinARunPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, LopatinASobelOperatorSEQ>(PPC_SETTINGS_lopatin_a_sobel_operator);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = LopatinARunPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunSobelPerfTests, LopatinARunPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace lopatin_a_sobel_operator
