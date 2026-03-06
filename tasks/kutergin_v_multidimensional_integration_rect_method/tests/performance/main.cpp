#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "../../common/include/common.hpp"
#include "../../omp/include/rect_method_omp.hpp"
#include "../../seq/include/rect_method_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kutergin_v_multidimensional_integration_rect_method {

class RectMethodPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  RectMethodPerfTests() = default;

 protected:
  void SetUp() override {
    input_data_.limits.resize(5, {0.0, 1.0});
    input_data_.n_steps = {40, 40, 40, 40, 40};
    input_data_.func = [](const std::vector<double> &x) {
      double sum = 0.0;
      for (double val : x) {
        sum += std::sin(val);
      }
      return sum;
    };
  }

  bool CheckTestOutputData([[maybe_unused]] OutType &output_data) final {
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
};

namespace {
TEST_P(RectMethodPerfTests, PerfTest) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, RectMethodSequential, RectMethodOMP>(
    PPC_SETTINGS_kutergin_v_multidimensional_integration_rect_method);

INSTANTIATE_TEST_SUITE_P(MultidimensionalIntegrationPerf, RectMethodPerfTests,
                         ppc::util::TupleToGTestValues(kAllPerfTasks), RectMethodPerfTests::CustomPerfTestName);

}  // namespace

}  // namespace kutergin_v_multidimensional_integration_rect_method
