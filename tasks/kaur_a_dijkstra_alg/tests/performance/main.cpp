#include <gtest/gtest.h>

#include <tuple>

#include "../../all/include/all.hpp"
#include "../../common/include/common.hpp"
#include "../../omp/include/omp.hpp"
#include "../../seq/include/seq.hpp"
#include "../../stl/include/stl.hpp"
#include "../../tbb/include/tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace kaur_a_dijkstra_alg {

class DijkstraPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_{};
  OutType expected_output_{};

  void SetUp() override {
    input_data_ = 50;
    expected_output_ = input_data_ * (input_data_ - 1) / 2;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return expected_output_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(DijkstraPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kSeqPerfTasks = ppc::util::MakeAllPerfTasks<InType, KaurADijkstraAlgSEQ>(PPC_SETTINGS_kaur_a_dijkstra_alg);
const auto kOmpPerfTasks = ppc::util::MakeAllPerfTasks<InType, KaurADijkstraAlgOMP>(PPC_SETTINGS_kaur_a_dijkstra_alg);
const auto kStlPerfTasks = ppc::util::MakeAllPerfTasks<InType, KaurADijkstraAlgSTL>(PPC_SETTINGS_kaur_a_dijkstra_alg);
const auto kTbbPerfTasks = ppc::util::MakeAllPerfTasks<InType, KaurADijkstraAlgTBB>(PPC_SETTINGS_kaur_a_dijkstra_alg);
const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, KaurADijkstraAlgALL>(PPC_SETTINGS_kaur_a_dijkstra_alg);

const auto kTestPerfTasksLists =
    std::tuple_cat(kSeqPerfTasks, kOmpPerfTasks, kStlPerfTasks, kTbbPerfTasks, kAllPerfTasks);
const auto kGtestValues = ppc::util::TupleToGTestValues(kTestPerfTasksLists);
const auto kPerfTestName = DijkstraPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(DijkstraSeqPerf, DijkstraPerfTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace kaur_a_dijkstra_alg
