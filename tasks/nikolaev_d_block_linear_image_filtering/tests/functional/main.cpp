#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "nikolaev_d_block_linear_image_filtering/common/include/common.hpp"
#include "nikolaev_d_block_linear_image_filtering/omp/include/ops_omp.hpp"
#include "nikolaev_d_block_linear_image_filtering/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace nikolaev_d_block_linear_image_filtering {

class NikolaevDBlockLinearImageFilteringFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_result_ = std::get<1>(params);
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

TEST_P(NikolaevDBlockLinearImageFilteringFuncTests, GaussBlur) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParam = {
    std::make_tuple(std::make_tuple(1, 1, std::vector<uint8_t>{0, 0, 0}), std::vector<uint8_t>{0, 0, 0}, "1pixel"),
    std::make_tuple(std::make_tuple(2, 1, std::vector<uint8_t>{200, 0, 0, 0, 200, 0}),
                    std::vector<uint8_t>{150, 50, 0, 50, 150, 0}, "2x1"),
    std::make_tuple(std::make_tuple(1, 2, std::vector<uint8_t>{255, 255, 255, 0, 0, 0}),
                    std::vector<uint8_t>{191, 191, 191, 64, 64, 64}, "1x2"),
    std::make_tuple(std::make_tuple(3, 3, std::vector<uint8_t>{0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 160, 160,
                                                               160, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}),
                    std::vector<uint8_t>{10, 10, 10, 20, 20, 20, 10, 10, 10, 20, 20, 20, 40, 40,
                                         40, 20, 20, 20, 10, 10, 10, 20, 20, 20, 10, 10, 10},
                    "3x3_center_point"),
    std::make_tuple(std::make_tuple(4, 1, std::vector<uint8_t>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}),
                    std::vector<uint8_t>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "4x1")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<NikolaevDBlockLinearImageFilteringSEQ, InType>(
                                               kTestParam, PPC_SETTINGS_nikolaev_d_block_linear_image_filtering),
                                           ppc::util::AddFuncTask<NikolaevDBlockLinearImageFilteringOMP, InType>(
                                               kTestParam, PPC_SETTINGS_nikolaev_d_block_linear_image_filtering));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    NikolaevDBlockLinearImageFilteringFuncTests::PrintFuncTestName<NikolaevDBlockLinearImageFilteringFuncTests>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, NikolaevDBlockLinearImageFilteringFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace nikolaev_d_block_linear_image_filtering
