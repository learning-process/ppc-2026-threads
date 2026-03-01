#include <gtest/gtest.h>

#include <chrono>
#include <cstddef>
#include <random>
#include <vector>

#include "fedoseev_linear_image_filtering_vertical/common/include/common.hpp"
#include "fedoseev_linear_image_filtering_vertical/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace fedoseev_linear_image_filtering_vertical {

class FedoseevPerfTest : public ppc::util::BaseRunPerfTests<Image, Image> {
 protected:
  void SetUp() override {
    const int size = 1024;
    input_.width = size;
    input_.height = size;
    input_.data.resize(static_cast<size_t>(size) * static_cast<size_t>(size));

    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed);
    std::uniform_int_distribution<int> dist(0, 255);
    for (auto &v : input_.data) {
      v = dist(gen);
    }
  }

  bool CheckTestOutputData(Image &output_data) override {
    return (output_data.width == input_.width && output_data.height == input_.height);
  }

  Image GetTestInputData() override {
    return input_;
  }

 private:
  Image input_;
};

TEST_P(FedoseevPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<Image, LinearImageFilteringVerticalSeq>(
    PPC_SETTINGS_fedoseev_linear_image_filtering_vertical);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = FedoseevPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerfTests, FedoseevPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace fedoseev_linear_image_filtering_vertical
