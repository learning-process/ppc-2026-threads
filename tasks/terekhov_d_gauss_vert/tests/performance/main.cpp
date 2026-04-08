#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <random>
#include <tuple>
#include <vector>

#include "terekhov_d_gauss_vert/common/include/common.hpp"
#include "terekhov_d_gauss_vert/omp/include/ops_omp.hpp"
#include "terekhov_d_gauss_vert/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace terekhov_d_gauss_vert {

class TerekhovDiGaussVertPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    const size_t total_pixels = 10000000;
    int img_size = static_cast<int>(std::sqrt(static_cast<double>(total_pixels)));
    input_data_.width = img_size;
    input_data_.height = img_size;
    input_data_.data.resize(static_cast<size_t>(input_data_.width) * static_cast<size_t>(input_data_.height));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 255);
    for (auto &pixel : input_data_.data) {
      pixel = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.width == input_data_.width && output_data.height == input_data_.height &&
           output_data.data.size() == input_data_.data.size();
  }

  InType GetTestInputData() final {
    return input_data_;
  }
  // 123344
 private:
  InType input_data_;
};

TEST_P(TerekhovDiGaussVertPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    std::tuple_cat(ppc::util::MakeAllPerfTasks<InType, TerekhovDGaussVertSEQ>(PPC_SETTINGS_terekhov_d_gauss_vert),
                   ppc::util::MakeAllPerfTasks<InType, TerekhovDGaussVertOMP>(PPC_SETTINGS_terekhov_d_gauss_vert));

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = TerekhovDiGaussVertPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerfTests, TerekhovDiGaussVertPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace terekhov_d_gauss_vert
