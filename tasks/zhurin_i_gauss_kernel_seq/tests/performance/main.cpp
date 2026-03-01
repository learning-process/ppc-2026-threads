#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <tuple>
#include <vector>

#include "util/include/perf_test_util.hpp"
#include "zhurin_i_gauss_kernel_seq/common/include/common.hpp"
#include "zhurin_i_gauss_kernel_seq/seq/include/ops_seq.hpp"

namespace zhurin_i_gauss_kernel_seq {

class ZhurinIGaussKernelPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    const int width = 512;
    const int height = 512;
    const int parts = 4;

    std::vector<std::vector<int>> img(height, std::vector<int>(width, 128));
    input_data_ = std::make_tuple(width, height, parts, img);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const auto &in = GetTestInputData();
    int expected_height = std::get<1>(in);
    int expected_width = std::get<0>(in);

    if (output_data.size() != static_cast<size_t>(expected_height)) {
      return false;
    }
    return std::ranges::all_of(
        output_data, [expected_width](const auto &row) { return row.size() == static_cast<size_t>(expected_width); });
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(ZhurinIGaussKernelPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ZhurinIGaussKernelSEQ>(PPC_SETTINGS_zhurin_i_gauss_kernel_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ZhurinIGaussKernelPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ZhurinIGaussKernelPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace zhurin_i_gauss_kernel_seq
