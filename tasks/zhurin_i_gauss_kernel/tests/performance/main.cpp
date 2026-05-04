#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <vector>

#include "util/include/perf_test_util.hpp"
#include "zhurin_i_gauss_kernel/common/include/common.hpp"
#include "zhurin_i_gauss_kernel/tbb/include/ops_tbb.hpp"

namespace zhurin_i_gauss_kernel {

class ZhurinIGaussKernelPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    const int width = 4096;
    const int height = 2048;
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

const auto kStlPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ZhurinIGaussKernelSTL>(PPC_SETTINGS_zhurin_i_gauss_kernel);
const auto kAllPerfTasks = std::tuple_cat(kStlPerfTasks);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ZhurinIGaussKernelPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ZhurinIGaussKernelPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace zhurin_i_gauss_kernel
