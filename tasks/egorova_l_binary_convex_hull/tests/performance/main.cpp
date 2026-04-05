#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

#include "egorova_l_binary_convex_hull/common/include/common.hpp"
#include "egorova_l_binary_convex_hull/omp/include/ops_omp.hpp"
#include "egorova_l_binary_convex_hull/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace egorova_l_binary_convex_hull {

class EgorovaLPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  static constexpr int kWidth = 1500;
  static constexpr int kHeight = 1500;
  InType input_data;
  int expected_components{};

  void SetUp() override {
    input_data.width = kWidth;
    input_data.height = kHeight;
    input_data.data.assign(static_cast<size_t>(kWidth) * static_cast<size_t>(kHeight), 0);

    expected_components = 0;
    for (int row = 100; row < kHeight - 100; row += 100) {
      for (int col = 100; col < kWidth - 100; col += 100) {
        for (int dy = 0; dy < 20; ++dy) {
          for (int dx = 0; dx < 20; ++dx) {
            const size_t index =
                (static_cast<size_t>(row + dy) * static_cast<size_t>(kWidth)) + static_cast<size_t>(col + dx);
            input_data.data[index] = 255;
          }
        }
        ++expected_components;
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.size() == static_cast<size_t>(expected_components);
  }

  InType GetTestInputData() final {
    return input_data;
  }
};

TEST_P(EgorovaLPerfTest, RunPerfModesSEQ) {
  ExecuteTest(GetParam());
}

namespace {
const auto kAllPerfTasksSEQ =
    ppc::util::MakeAllPerfTasks<InType, BinaryConvexHullSEQ>(PPC_SETTINGS_egorova_l_binary_convex_hull);
INSTANTIATE_TEST_SUITE_P(RunModeTestsSEQ, EgorovaLPerfTest, ppc::util::TupleToGTestValues(kAllPerfTasksSEQ),
                         EgorovaLPerfTest::CustomPerfTestName);
}  // namespace

TEST_P(EgorovaLPerfTest, RunPerfModesOMP) {
  ExecuteTest(GetParam());
}

namespace {
const auto kAllPerfTasksOMP =
    ppc::util::MakeAllPerfTasks<InType, BinaryConvexHullOMP>(PPC_SETTINGS_egorova_l_binary_convex_hull);
INSTANTIATE_TEST_SUITE_P(RunModeTestsOMP, EgorovaLPerfTest, ppc::util::TupleToGTestValues(kAllPerfTasksOMP),
                         EgorovaLPerfTest::CustomPerfTestName);
}  // namespace

}  // namespace egorova_l_binary_convex_hull
