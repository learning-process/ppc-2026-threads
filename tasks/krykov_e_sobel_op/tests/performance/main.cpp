#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>

#include "krykov_e_sobel_op/common/include/common.hpp"
#include "krykov_e_sobel_op/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace krykov_e_sobel_op {

class KrykovERunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  void SetUp() override {
    const int size = 512;

    Image img;
    img.width = size;
    img.height = size;
    img.data.resize(static_cast<size_t>(size) * static_cast<size_t>(size));

    for (int row = 0; row < size; ++row) {
      for (int col = 0; col < size; ++col) {
        auto value = static_cast<uint8_t>((255 * col) / (size - 1));
        img.data[(row * size) + col] = {.r = value, .g = value, .b = value};
      }
    }

    input_data_ = img;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return static_cast<int>(output_data.size()) == input_data_.width * input_data_.height;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(KrykovERunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, KrykovESobelOpSEQ>(PPC_SETTINGS_krykov_e_sobel_op);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KrykovERunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KrykovERunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace krykov_e_sobel_op
