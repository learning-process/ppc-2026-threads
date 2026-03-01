#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ilin_a_algorithm_graham/common/include/common.hpp"
#include "ilin_a_algorithm_graham/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace ilin_a_algorithm_graham {

class IlinAGrahamFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return "id_" + std::to_string(std::get<0>(test_param));
  }

 protected:
  void SetUp() override {
    auto params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<1>(params);
    expected_ = std::get<2>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_.size()) {
      return false;
    }

    for (size_t i = 0; i < output_data.size(); ++i) {
      if (std::abs(output_data[i].x - expected_[i].x) > 1e-6 || std::abs(output_data[i].y - expected_[i].y) > 1e-6) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_;
};

namespace {

InputData MakeInput(std::vector<Point> points) {
  return InputData{.points = std::move(points)};
}

const std::array<TestType, 3> kTestCases = {
    std::make_tuple(1, MakeInput({{0, 0}, {1, 0}, {0, 1}, {1, 1}}), OutType{{0, 0}, {1, 0}, {1, 1}, {0, 1}}),

    std::make_tuple(2, MakeInput({{0, 0}, {2, 0}, {1, 1}, {0, 2}, {2, 2}, {1, 3}}),
                    OutType{{0, 0}, {2, 0}, {2, 2}, {1, 3}, {0, 2}}),

    std::make_tuple(3, MakeInput({{0, 0}, {3, 0}, {1, 1}, {2, 1}, {0, 3}, {3, 3}, {1, 2}, {2, 2}}),
                    OutType{{0, 0}, {3, 0}, {3, 3}, {0, 3}})};

const auto kTestTasksList =
    ppc::util::AddFuncTask<IlinAGrahamSEQ, InType>(kTestCases, PPC_SETTINGS_ilin_a_algorithm_graham);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = IlinAGrahamFuncTests::PrintFuncTestName<IlinAGrahamFuncTests>;

INSTANTIATE_TEST_SUITE_P(GrahamTests, IlinAGrahamFuncTests, kGtestValues, kTestName);

TEST_P(IlinAGrahamFuncTests, BuildConvexHull) {
  ExecuteTest(GetParam());
}

}  // namespace

}  // namespace ilin_a_algorithm_graham
