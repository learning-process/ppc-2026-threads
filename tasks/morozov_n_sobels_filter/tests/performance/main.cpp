#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <random>

#include "morozov_n_sobels_filter/common/include/common.hpp"
#include "morozov_n_sobels_filter/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace morozov_n_sobels_filter {

class MorozovNRunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;
  OutType correct_image_;

  void SetUp() override {
    std::size_t height = 5000;
    std::size_t width = 5000;
    int seed = 777;

    GenerateImage(height, width, seed);
    correct_image_.pixels.resize(height * width);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.pixels.size() == correct_image_.pixels.size();
  }

  void GenerateImage(size_t height, size_t width, int seed) {
    input_data_.height = height;
    input_data_.width = width;
    input_data_.pixels.resize(height * width, 0);

    std::mt19937 gen(seed);
    std::uniform_int_distribution pixel_dis(0, 255);

    for (unsigned char &pixel : input_data_.pixels) {
      pixel = static_cast<std::uint8_t>(pixel_dis(gen));
    }
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(MorozovNRunPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, MorozovNSobelsFilterSEQ>(PPC_SETTINGS_morozov_n_sobels_filter);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = MorozovNRunPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunSobelPerfTests, MorozovNRunPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace morozov_n_sobels_filter
