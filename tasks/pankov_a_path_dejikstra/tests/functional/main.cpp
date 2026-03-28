#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <tuple>

#include "pankov_a_path_dejikstra/common/include/common.hpp"
#include "pankov_a_path_dejikstra/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace pankov_a_path_dejikstra {

class PankovAPathDejikstraRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    const TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    const int case_id = std::get<0>(params);

    if (case_id == 0) {
      input_data_ = InType{5, 0, {{0, 1, 4}, {0, 2, 1}, {2, 1, 2}, {1, 3, 1}, {2, 3, 5}, {3, 4, 3}}};
      expected_output_ = OutType{0, 3, 1, 4, 7};
      return;
    }
    if (case_id == 1) {
      input_data_ = InType{6, 0, {{0, 1, 7}, {0, 2, 9}, {1, 3, 10}, {2, 3, 1}, {3, 4, 2}}};
      expected_output_ = OutType{0, 7, 9, 10, 12, kInfinity};
      return;
    }
    if (case_id == 2) {
      input_data_ = InType{5, 2, {{2, 0, 5}, {2, 1, 2}, {1, 3, 4}, {0, 3, 1}, {3, 4, 3}}};
      expected_output_ = OutType{5, 2, 0, 6, 9};
      return;
    }
    throw std::runtime_error("Unknown functional case id");
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == expected_output_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
  OutType expected_output_;
};

namespace {

TEST_P(PankovAPathDejikstraRunFuncTestsThreads, DijkstraShortestPathFromSource) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(0, "basic_graph"), std::make_tuple(1, "unreachable_vertex"),
                                            std::make_tuple(2, "non_zero_source")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<PankovAPathDejikstraSEQ, InType>(kTestParam, PPC_SETTINGS_pankov_a_path_dejikstra);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    PankovAPathDejikstraRunFuncTestsThreads::PrintFuncTestName<PankovAPathDejikstraRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, PankovAPathDejikstraRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace pankov_a_path_dejikstra
