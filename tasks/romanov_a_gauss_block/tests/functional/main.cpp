#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "romanov_a_gauss_block/common/include/common.hpp"
#include "romanov_a_gauss_block/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace romanov_a_gauss_block {

class RomanovARunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<4>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::make_tuple(std::get<0>(params), std::get<1>(params), std::get<2>(params));
    expected_result_ = std::get<3>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (expected_result_ == output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_result_;
};

namespace {

TEST_P(RomanovARunFuncTestsThreads, GaussFilter) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 7> kTestParam = {
    std::make_tuple(1, 1, std::vector<uint8_t>{0, 0, 0}, std::vector<uint8_t>{0, 0, 0}, "Test1"),
    std::make_tuple(2, 1, std::vector<uint8_t>{224, 143, 37, 137, 16, 22}, std::vector<uint8_t>{73, 38, 12, 62, 22, 10},
                    "Test2"),
    std::make_tuple(1, 2, std::vector<uint8_t>{224, 143, 37, 137, 16, 22}, std::vector<uint8_t>{73, 38, 12, 62, 22, 10},
                    "Test3"),
    std::make_tuple(2, 2, std::vector<uint8_t>{100, 150, 200, 50, 60, 70, 80, 90, 100, 120, 130, 140},
                    std::vector<uint8_t>{49, 64, 80, 45, 56, 66, 51, 61, 72, 53, 61, 69}, "Test4"),
    std::make_tuple(
        3, 3, std::vector<uint8_t>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 201, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        std::vector<uint8_t>{13, 0, 0, 25, 0, 0, 13, 0, 0, 25, 0, 0, 50, 0, 0, 25, 0, 0, 13, 0, 0, 25, 0, 0, 13, 0, 0},
        "Test5"),
    std::make_tuple(2, 3, std::vector<uint8_t>{10, 10, 10, 20, 20, 20, 30, 30, 30, 40, 40, 40, 50, 50, 50, 60, 60, 60},
                    std::vector<uint8_t>{11, 11, 11, 13, 13, 13, 25, 25, 25, 28, 28, 28, 26, 26, 26, 28, 28, 28},
                    "Test6"),
    std::make_tuple(4, 1, std::vector<uint8_t>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    std::vector<uint8_t>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "Test7")};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<RomanovAGaussBlockSEQ, InType>(kTestParam, PPC_SETTINGS_example_threads));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = RomanovARunFuncTestsThreads::PrintFuncTestName<RomanovARunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(GaussFilter, RomanovARunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace romanov_a_gauss_block
