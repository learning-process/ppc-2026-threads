#include <gtest/gtest.h>
#include <mpi.h>

#include <cmath>
#include <tuple>
#include <vector>

#include "shkrebko_m_calc_of_integral_rect/all/include/ops_all.hpp"
#include "shkrebko_m_calc_of_integral_rect/common/include/common.hpp"
#include "shkrebko_m_calc_of_integral_rect/omp/include/ops_omp.hpp"
#include "shkrebko_m_calc_of_integral_rect/seq/include/ops_seq.hpp"
#include "shkrebko_m_calc_of_integral_rect/stl/include/ops_stl.hpp"
#include "shkrebko_m_calc_of_integral_rect/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"
#include "util/include/util.hpp"

namespace shkrebko_m_calc_of_integral_rect {

class ShkrebkoMRunPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  static constexpr int kN = 200;

 protected:
  void SetUp() override {
    input_data_.func = [](const std::vector<double> &x) { return std::sin(x[0]) + std::cos(x[1]) + x[2]; };
    input_data_.limits = {{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}};
    input_data_.n_steps = {kN, kN, kN};

    int rank = 0;
    if (ppc::util::IsUnderMpirun()) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
    if (rank == 0) {
      expected_ = (1.0 - std::cos(1.0)) + std::sin(1.0) + 0.5;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    int rank = 0;
    if (ppc::util::IsUnderMpirun()) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
    if (rank != 0) {
      return true;
    }
    const double eps = 1e-3;
    return std::fabs(output_data - expected_) <= eps;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_ = 0.0;
};

TEST_P(ShkrebkoMRunPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {
const auto kAllPerfTasks = std::tuple_cat(
    ppc::util::MakeAllPerfTasks<InType, ShkrebkoMCalcOfIntegralRectSEQ>(PPC_SETTINGS_shkrebko_m_calc_of_integral_rect),
    ppc::util::MakeAllPerfTasks<InType, ShkrebkoMCalcOfIntegralRectOMP>(PPC_SETTINGS_shkrebko_m_calc_of_integral_rect),
    ppc::util::MakeAllPerfTasks<InType, ShkrebkoMCalcOfIntegralRectTBB>(PPC_SETTINGS_shkrebko_m_calc_of_integral_rect),
    ppc::util::MakeAllPerfTasks<InType, ShkrebkoMCalcOfIntegralRectSTL>(PPC_SETTINGS_shkrebko_m_calc_of_integral_rect),
    ppc::util::MakeAllPerfTasks<InType, ShkrebkoMCalcOfIntegralRectALL>(PPC_SETTINGS_shkrebko_m_calc_of_integral_rect));

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ShkrebkoMRunPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ShkrebkoMRunPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace shkrebko_m_calc_of_integral_rect
