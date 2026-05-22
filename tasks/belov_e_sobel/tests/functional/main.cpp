#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <cstdint>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "belov_e_sobel/all/include/ops_all.hpp"
#include "belov_e_sobel/common/include/common.hpp"
#include "belov_e_sobel/omp/include/ops_omp.hpp"
#include "belov_e_sobel/seq/include/ops_seq.hpp"
#include "belov_e_sobel/stl/include/ops_stl.hpp"
#include "belov_e_sobel/tbb/include/ops_tbb.hpp"
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
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam()) + ".jpg";
    std::string abs_path = ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_belov_e_sobel), params);

    int w;
    int h;
    int c;
    uint8_t *raw_data = stbi_load(abs_path.c_str(), &w, &h, &c, 1);
    if (raw_data == nullptr) {
      std::cout << abs_path << std::endl;
    }

    std::vector<uint8_t> input_vec(raw_data, raw_data + (static_cast<ptrdiff_t>(w * h)));
    stbi_image_free(raw_data);

    input_data_ = std::make_tuple(input_vec, w, h);
  }
  InType GetTestInputData() final {
    return input_data_;
  }
  bool CheckTestOutputData(OutType &output_data) final {
    std::string params =
        "ref" + std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam()) + ".jpg";
    std::string abs_path = ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_belov_e_sobel), params);

    int w;
    int h;
    int c;
    uint8_t *raw_data = stbi_load(abs_path.c_str(), &w, &h, &c, 1);

    std::vector<uint8_t> ref_vector(raw_data, raw_data + (static_cast<ptrdiff_t>(w * h)));
    stbi_image_free(raw_data);

    auto [out_vector, out_w, out_h] = output_data;
    if (w != out_w || h != out_h || ref_vector.size() != out_vector.size()) {
      return false;
    }

    const int kPixelThreshold = 1;

    size_t structural_errors = 0;

    for (size_t i = 0; i < ref_vector.size(); ++i) {
      int diff = std::abs(static_cast<int>(ref_vector[i]) - static_cast<int>(out_vector[i]));
      if (diff > kPixelThreshold) {
        structural_errors++;
      }
    }

    size_t max_allowed_errors = ref_vector.size() * 0.0005;

    return structural_errors <= max_allowed_errors;
  }
};

namespace {

TEST_P(BelovESobelFuncTests, Sobel) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 2> kTestParam = {"test1", "test2"};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<BelovESobelALL, InType>(kTestParam, PPC_SETTINGS_belov_e_sobel),
                   ppc::util::AddFuncTask<BelovESobelOMP, InType>(kTestParam, PPC_SETTINGS_belov_e_sobel),
                   ppc::util::AddFuncTask<BelovESobelSEQ, InType>(kTestParam, PPC_SETTINGS_belov_e_sobel),
                   ppc::util::AddFuncTask<BelovESobelSTL, InType>(kTestParam, PPC_SETTINGS_belov_e_sobel),
                   ppc::util::AddFuncTask<BelovESobelTBB, InType>(kTestParam, PPC_SETTINGS_belov_e_sobel));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = BelovESobelFuncTests::PrintFuncTestName<BelovESobelFuncTests>;

INSTANTIATE_TEST_SUITE_P(SobelTests, BelovESobelFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace belov_e_sobel
