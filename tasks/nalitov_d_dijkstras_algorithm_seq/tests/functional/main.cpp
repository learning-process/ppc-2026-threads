#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

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
    double tol = 5.0 / std::sqrt(static_cast<double>(input_data_.num_samples));
    return std::abs(output_data - ExactValue(input_data_)) <= tol;
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

const std::array<TestType, 3> kTestParam = {std::make_tuple(5, "5"), std::make_tuple(7, "7"), std::make_tuple(9, "9")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<NalitovDDijkstrasAlgorithmSeq, InType>(kTestParam, PPC_SETTINGS_nalitov_d_dijkstras_algorithm_seq);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = NalitovDDijkstrasAlgorithmSeqFuncTests::PrintFuncTestName<NalitovDDijkstrasAlgorithmSeqFuncTests>;

INSTANTIATE_TEST_SUITE_P(DijkstraAlgorithmTests, NalitovDDijkstrasAlgorithmSeqFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace nalitov_d_dijkstras_algorithm_seq