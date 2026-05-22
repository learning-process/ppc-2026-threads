#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "../../all/include/all.hpp"
#include "../../common/include/common.hpp"
#include "../../omp/include/omp.hpp"
#include "../../seq/include/seq.hpp"
#include "../../stl/include/stl.hpp"
#include "../../tbb/include/tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kaur_a_dijkstra_alg {

class KaurARunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_output_ = input_data_ * (input_data_ - 1) / 2;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return expected_output_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = 0;
  OutType expected_output_ = 0;
};

namespace {

TEST_P(KaurARunFuncTestsThreads, DijkstraFromParams) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(3, "3"), std::make_tuple(5, "5"), std::make_tuple(7, "7")};

const auto kSeqTasks = ppc::util::AddFuncTask<KaurADijkstraAlgSEQ, InType>(kTestParam, PPC_SETTINGS_kaur_a_dijkstra_alg);
const auto kOmpTasks = ppc::util::AddFuncTask<KaurADijkstraAlgOMP, InType>(kTestParam, PPC_SETTINGS_kaur_a_dijkstra_alg);
const auto kStlTasks = ppc::util::AddFuncTask<KaurADijkstraAlgSTL, InType>(kTestParam, PPC_SETTINGS_kaur_a_dijkstra_alg);
const auto kTbbTasks = ppc::util::AddFuncTask<KaurADijkstraAlgTBB, InType>(kTestParam, PPC_SETTINGS_kaur_a_dijkstra_alg);
const auto kAllTasks = ppc::util::AddFuncTask<KaurADijkstraAlgALL, InType>(kTestParam, PPC_SETTINGS_kaur_a_dijkstra_alg);
const auto kTestTasksList = std::tuple_cat(kSeqTasks, kOmpTasks, kStlTasks, kTbbTasks, kAllTasks);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kPerfTestName = KaurARunFuncTestsThreads::PrintFuncTestName<KaurARunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(DijkstraSeqTests, KaurARunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kaur_a_dijkstra_alg
