#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "nalitov_d_dijkstras_algorithm_seq/common/include/common.hpp"
#include "nalitov_d_dijkstras_algorithm_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace nalitov_d_dijkstras_algorithm_seq {

class NalitovDDijkstrasAlgorithmSeqFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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

TEST_P(NalitovDDijkstrasAlgorithmSeqFuncTests, AlgorithmIntegration) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(2, "2"), std::make_tuple(4, "4"), std::make_tuple(6, "6")};

const auto kTestTasksList = ppc::util::AddFuncTask<NalitovDDijkstrasAlgorithmSeq, InType>(
    kTestParam, PPC_SETTINGS_nalitov_d_dijkstras_algorithm_seq);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    NalitovDDijkstrasAlgorithmSeqFuncTests::PrintFuncTestName<NalitovDDijkstrasAlgorithmSeqFuncTests>;

INSTANTIATE_TEST_SUITE_P(DijkstraAlgorithmTests, NalitovDDijkstrasAlgorithmSeqFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace nalitov_d_dijkstras_algorithm_seq
