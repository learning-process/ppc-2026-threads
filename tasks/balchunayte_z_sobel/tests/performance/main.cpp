#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>

#include "balchunayte_z_sobel/common/include/common.hpp"
#include "balchunayte_z_sobel/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace balchunayte_z_sobel {

class BalchunayteZRunPerfTestSEQ : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    const int image_size = 512;

    Image input_image;
    input_image.width = image_size;
    input_image.height = image_size;
    input_image.data.resize(static_cast<size_t>(image_size) * static_cast<size_t>(image_size));

    for (int row_index = 0; row_index < image_size; ++row_index) {
      for (int col_index = 0; col_index < image_size; ++col_index) {
        const auto intensity_value = static_cast<uint8_t>((255 * col_index) / (image_size - 1));
        input_image.data[static_cast<size_t>((row_index * image_size) + col_index)] =
            Pixel{intensity_value, intensity_value, intensity_value};
      }
    }

    input_data_ = input_image;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return static_cast<int>(output_data.size()) == input_data_.width * input_data_.height;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
};

TEST_P(BalchunayteZRunPerfTestSEQ, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, BalchunayteZSobelOpSEQ>(PPC_SETTINGS_balchunayte_z_sobel);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = BalchunayteZRunPerfTestSEQ::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, BalchunayteZRunPerfTestSEQ, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace balchunayte_z_sobel
