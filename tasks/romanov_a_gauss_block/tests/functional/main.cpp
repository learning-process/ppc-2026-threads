#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
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
    return (std::get<2>(input_data_) == output_data);
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

const std::array<TestType, 1> kTestParam = {
  std::make_tuple(1, 1, std::vector<uint8_t>{0, 0, 0}, std::vector<uint8_t>{0, 0, 0}, "TestNo1")
};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<RomanovAGaussBlockSEQ, InType>(kTestParam, PPC_SETTINGS_example_threads));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = RomanovARunFuncTestsThreads::PrintFuncTestName<RomanovARunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(GaussFilter, RomanovARunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace romanov_a_gauss_block
