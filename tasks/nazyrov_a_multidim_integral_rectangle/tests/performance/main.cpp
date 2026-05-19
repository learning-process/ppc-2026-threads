#include <gtest/gtest.h>
#include <mpi.h>

#include <cmath>
#include <tuple>
#include <vector>

#include "nazyrov_a_multidim_integral_rectangle/common/include/common.hpp"
#include "nazyrov_a_multidim_integral_rectangle/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"
#include "util/include/util.hpp"

namespace nazyrov_a_multidim_integral_rectangle {

class NazyrovAMultidimIntegralRectanglePerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kDim = 3;
  static constexpr int kSteps = 300;
  InType input_data_;

 public:
  void SetUp() override {
    Func func = [](const std::vector<double> &x) {
      double s = 0.0;
      for (double xi : x) {
        s += xi * xi;
      }
      return s;
    };
    input_data_ = std::make_tuple(func, Bounds(kDim, {0.0, 1.0}), kSteps);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    int rank = 0;
    if (ppc::util::IsUnderMpirun()) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
    if (rank != 0) {
      return true;
    }
    // f = sum(x_i^2) over [0,1]^3, expected = 1.0
    return std::abs(output_data - 1.0) < 1e-3;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(NazyrovAMultidimIntegralRectanglePerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, NazyrovAMultidimIntegralRectangleSeq>(
    PPC_SETTINGS_nazyrov_a_multidim_integral_rectangle);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = NazyrovAMultidimIntegralRectanglePerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, NazyrovAMultidimIntegralRectanglePerfTest, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace nazyrov_a_multidim_integral_rectangle
