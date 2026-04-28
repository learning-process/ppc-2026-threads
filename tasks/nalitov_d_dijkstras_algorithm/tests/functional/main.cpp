#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "nalitov_d_dijkstras_algorithm/common/include/common.hpp"
#include "nalitov_d_dijkstras_algorithm/omp/include/ops_omp.hpp"
#include "nalitov_d_dijkstras_algorithm/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace nalitov_d_dijkstras_algorithm {

class NalitovDDijkstrasAlgorithmFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    OutType expected_output = input_data_ * (input_data_ - 1) / 2;
    return expected_output == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
};

namespace {

TEST_P(NalitovDDijkstrasAlgorithmFuncTests, AlgorithmIntegration) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(2, "2"), std::make_tuple(4, "4"), std::make_tuple(6, "6")};

const auto kTestTasksListSeq = ppc::util::AddFuncTask<NalitovDDijkstrasAlgorithmSeq, InType>(
    kTestParam, PPC_SETTINGS_nalitov_d_dijkstras_algorithm);
const auto kTestTasksListOmp = ppc::util::AddFuncTask<NalitovDDijkstrasAlgorithmOmp, InType>(
    kTestParam, PPC_SETTINGS_nalitov_d_dijkstras_algorithm);
const auto kTestTasksList = std::tuple_cat(kTestTasksListSeq, kTestTasksListOmp);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = NalitovDDijkstrasAlgorithmFuncTests::PrintFuncTestName<NalitovDDijkstrasAlgorithmFuncTests>;

INSTANTIATE_TEST_SUITE_P(DijkstraAlgorithmTests, NalitovDDijkstrasAlgorithmFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace nalitov_d_dijkstras_algorithm
