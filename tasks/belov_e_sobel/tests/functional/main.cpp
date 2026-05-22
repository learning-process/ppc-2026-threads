#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <cstdint>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "util/include/func_test_util.hpp"

namespace belov_e_sobel {

class BelovESobelFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 private:
  InType input_data_;

 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return test_param;
  }

 protected:
  static void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    std::string abs_path = ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_belov_e_sobel), params);

    int w;
    int h;
    int c;
    uint8_t *raw_data = stbi_load(abs_path.c_str(), &w, &h, &c, 1);

    std::vector<uint8_t> input_vec(raw_data, raw_data + (static_cast<ptrdiff_t>(w * h)));
    stbi_image_free(raw_data);

    input_data_ = std::make_tuple(input_vec, w, h);
  }
  InType GetTestInputData() final {
    return input_data_;
  }
  static bool CheckTestOutputData(OutType &output_data) final {
    std::string params =
        "ref" + std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    std::string abs_path = ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_belov_e_sobel), params);

    int w;
    int h;
    int c;
    uint8_t *raw_data = stbi_load(abs_path.c_str(), &w, &h, &c, 1);

    std::vector<uint8_t> ref_vector(raw_data, raw_data + (static_cast<ptrdiff_t>(w * h)));
    stbi_image_free(raw_data);

    return (std::make_tuple(ref_vector, w, h) == output_data);
  }
};

namespace {

TEST_P(BelovESobelFuncTests, Sobel) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 2> kTestParam = {"test1.jpg", "test2.jpg"};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<BelovESobelALL, InType>(kTestParam, PPC_ID_belov_e_sobel),
                   ppc::util::AddFuncTask<BelovESobelOMP, InType>(kTestParam, PPC_ID_belov_e_sobel),
                   ppc::util::AddFuncTask<BelovESobelSEQ, InType>(kTestParam, PPC_ID_belov_e_sobel),
                   ppc::util::AddFuncTask<BelovESobelSTL, InType>(kTestParam, PPC_ID_belov_e_sobel),
                   ppc::util::AddFuncTask<BelovESobelTBB, InType>(kTestParam, PPC_ID_belov_e_sobel));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = BelovESobelFuncTests::PrintFuncTestName<BelovESobelFuncTests>;

INSTANTIATE_TEST_SUITE_P(SobelTests, BelovESobelFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace belov_e_sobel
