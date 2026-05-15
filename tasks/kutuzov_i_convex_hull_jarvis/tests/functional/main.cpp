#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "kutuzov_i_convex_hull_jarvis/common/include/common.hpp"
#include "kutuzov_i_convex_hull_jarvis/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kutuzov_i_convex_hull_jarvis {

class KutuzovIRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param).size()) + "_" + std::to_string(std::get<1>(test_param).size());
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_output_ = std::get<1>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (output_data == expected_output_);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_output_;
};

namespace {

TEST_P(KutuzovIRunFuncTestsThreads, MatmulFromPic) {
  ExecuteTest(GetParam());
}

// Square
const InType kTestInput1 = {{0, 0}, {4, 0}, {4, 4}, {0, 4}, {2, 2}, {1, 1}, {3, 3}, {2, 3}};
const OutType kTestOutput1 = {{0, 0}, {4, 0}, {4, 4}, {0, 4}};

// Line
const InType kTestInput2 = {{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5}};
const OutType kTestOutput2 = {{0, 0}, {5, 5}};

// Triangle
const InType kTestInput3 = {{0, 0}, {2, 2}, {3, 1}, {4, 2}, {3, 4}, {6, 0}};
const OutType kTestOutput3 = {{0, 0}, {6, 0}, {3, 4}};

const std::array<TestType, 3> kTestParam = {std::make_tuple(kTestInput1, kTestOutput1),
                                            std::make_tuple(kTestInput2, kTestOutput2),
                                            std::make_tuple(kTestInput3, kTestOutput3)};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<KutuzovITestConvexHullSEQ, InType>(kTestParam, PPC_SETTINGS_kutuzov_i_convex_hull_jarvis));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = KutuzovIRunFuncTestsThreads::PrintFuncTestName<KutuzovIRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, KutuzovIRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kutuzov_i_convex_hull_jarvis
