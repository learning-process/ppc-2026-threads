#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "kolotukhin_a_gaussian_blur/common/include/common.hpp"
#include "kolotukhin_a_gaussian_blur/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kolotukhin_a_gaussian_blur {

class KolotukhinAGaussinBlureFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    return get<1>(params) == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(KolotukhinAGaussinBlureFuncTests, GaussianBlure) {
  ExecuteTest(GetParam());
}

// [TEST CASE] pixel: corner
std::vector<std::uint8_t> test_1 = {255, 0, 0, 0, 0, 0, 0, 0, 0};
std::vector<std::uint8_t> expect_1 = {143, 47, 0, 47, 15, 0, 0, 0, 0};

// [TEST CASE] pixel: border
std::vector<std::uint8_t> test_2 = {0, 0, 0, 255, 0, 0, 0, 0, 0};
std::vector<std::uint8_t> expect_2 = {47, 15, 0, 95, 31, 0, 47, 15, 0};

// [TEST CASE] pixel: inside
std::vector<std::uint8_t> test_3 = {0, 0, 0, 0, 255, 0, 0, 0, 0};
std::vector<std::uint8_t> expect_3 = {15, 31, 15, 31, 63, 31, 15, 31, 15};

// [TEST CASE] pixels
std::vector<std::uint8_t> test_4 = {15, 15, 15, 0, 90, 87, 42, 1, 12, 13, 134, 12};
std::vector<std::uint8_t> expect_4 = {33, 30, 19, 5, 51, 52, 42, 17, 31, 51, 65, 34};

const std::array<TestType, 4> kTestParam = {std::make_tuple(InType{test_1, 3, 3}, expect_1, "test_corner"),
                                            std::make_tuple(InType{test_2, 3, 3}, expect_2, "test_border"),
                                            std::make_tuple(InType{test_3, 3, 3}, expect_3, "test_inside"),
                                            std::make_tuple(InType{test_4, 4, 3}, expect_4, "test_common")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<KolotukhinAGaussinBlureSEQ, InType>(kTestParam, PPC_SETTINGS_kolotukhin_a_gaussian_blur);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = KolotukhinAGaussinBlureFuncTests::PrintFuncTestName<KolotukhinAGaussinBlureFuncTests>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, KolotukhinAGaussinBlureFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kolotukhin_a_gaussian_blur
