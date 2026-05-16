#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "nalitov_d_dijkstras_algorithm/all/include/ops_all.hpp"
#include "nalitov_d_dijkstras_algorithm/common/include/common.hpp"
#include "nalitov_d_dijkstras_algorithm/omp/include/ops_omp.hpp"
#include "nalitov_d_dijkstras_algorithm/seq/include/ops_seq.hpp"
#include "nalitov_d_dijkstras_algorithm/stl/include/ops_stl.hpp"
#include "nalitov_d_dijkstras_algorithm/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace nalitov_d_dijkstras_algorithm {

namespace {

InType MakeStarFromZero(int n) {
  InType g;
  g.n = n;
  g.source = 0;
  g.arcs.reserve(static_cast<std::size_t>(std::max(0, n - 1)));
  for (int i = 1; i < n; ++i) {
    g.arcs.push_back(Arc{.from = 0, .to = i, .weight = i});
  }
  return g;
}

InType MakeTriangleShortcutGraph() {
  InType g;
  g.n = 3;
  g.source = 0;
  g.arcs = {Arc{.from = 0, .to = 1, .weight = 1}, Arc{.from = 1, .to = 2, .weight = 1},
            Arc{.from = 0, .to = 2, .weight = 5}};
  return g;
}

}  // namespace

class NalitovDDijkstrasAlgorithmFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_output_ = std::get<1>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return expected_output_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
  OutType expected_output_{};
};

namespace {

TEST_P(NalitovDDijkstrasAlgorithmFuncTests, AlgorithmIntegration) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 4> kTestParam = {
    std::make_tuple(MakeStarFromZero(2), 1, "star_n2"),
    std::make_tuple(MakeStarFromZero(4), 6, "star_n4"),
    std::make_tuple(MakeStarFromZero(6), 15, "star_n6"),
    std::make_tuple(MakeTriangleShortcutGraph(), 3, "triangle_shortcut"),
};

const auto kTestTasksListSeq = ppc::util::AddFuncTask<NalitovDDijkstrasAlgorithmSeq, InType>(
    kTestParam, PPC_SETTINGS_nalitov_d_dijkstras_algorithm);
const auto kTestTasksListOmp = ppc::util::AddFuncTask<NalitovDDijkstrasAlgorithmOmp, InType>(
    kTestParam, PPC_SETTINGS_nalitov_d_dijkstras_algorithm);
const auto kTestTasksListTbb = ppc::util::AddFuncTask<NalitovDDijkstrasAlgorithmTBB, InType>(
    kTestParam, PPC_SETTINGS_nalitov_d_dijkstras_algorithm);
const auto kTestTasksListStl = ppc::util::AddFuncTask<NalitovDDijkstrasAlgorithmSTL, InType>(
    kTestParam, PPC_SETTINGS_nalitov_d_dijkstras_algorithm);
const auto kTestTasksListAll = ppc::util::AddFuncTask<NalitovDDijkstrasAlgorithmALL, InType>(
    kTestParam, PPC_SETTINGS_nalitov_d_dijkstras_algorithm);

const auto kTestTasksList =
    std::tuple_cat(kTestTasksListSeq, kTestTasksListOmp, kTestTasksListTbb, kTestTasksListStl, kTestTasksListAll);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = NalitovDDijkstrasAlgorithmFuncTests::PrintFuncTestName<NalitovDDijkstrasAlgorithmFuncTests>;

INSTANTIATE_TEST_SUITE_P(DijkstraAlgorithmTests, NalitovDDijkstrasAlgorithmFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace nalitov_d_dijkstras_algorithm
