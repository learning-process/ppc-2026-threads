#include <gtest/gtest.h>
#include <mpi.h>

#include <cmath>

#include "makovskiy_i_graham_hull/all/include/ops_all.hpp"
#include "makovskiy_i_graham_hull/common/include/common.hpp"
#include "makovskiy_i_graham_hull/omp/include/ops_omp.hpp"
#include "makovskiy_i_graham_hull/seq/include/ops_seq.hpp"
#include "makovskiy_i_graham_hull/stl/include/ops_stl.hpp"
#include "makovskiy_i_graham_hull/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"
#include "util/include/util.hpp"

namespace makovskiy_i_graham_hull {

class MakovskiyIGrahamHullRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    const int num_points = 500000;
    input_data_.clear();
    input_data_.resize(num_points);
    for (int i = 0; i < num_points; ++i) {
      input_data_[i] = {.x = std::sin(static_cast<double>(i)) * 100.0, .y = std::cos(static_cast<double>(i)) * 100.0};
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (ppc::util::IsUnderMpirun()) {
      int rank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank != 0) {
        return true;
      }
    }
    return output_data.size() >= 3;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(MakovskiyIGrahamHullRunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ConvexHullGrahamSEQ, ConvexHullGrahamOMP, ConvexHullGrahamTBB,
                                ConvexHullGrahamSTL, ConvexHullGrahamALL>(PPC_SETTINGS_makovskiy_i_graham_hull);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = MakovskiyIGrahamHullRunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, MakovskiyIGrahamHullRunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace makovskiy_i_graham_hull
