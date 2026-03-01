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

#include "zaharov_g_linear_contrast_stretch/all/include/ops_all.hpp"
#include "zaharov_g_linear_contrast_stretch/common/include/common.hpp"
#include "zaharov_g_linear_contrast_stretch/omp/include/ops_omp.hpp"
#include "zaharov_g_linear_contrast_stretch/seq/include/ops_seq.hpp"
#include "zaharov_g_linear_contrast_stretch/stl/include/ops_stl.hpp"
#include "zaharov_g_linear_contrast_stretch/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace zaharov_g_linear_contrast_stretch {

class ZaharovGRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    int width = -1;
    int height = -1;
    int channels = -1;
    std::vector<uint8_t> img;
    // Read image in RGB to ensure consistent channel count
    {
      std::string abs_path = ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_zaharov_g_linear_contrast_stretch), "pic.ppm");
      auto *data = stbi_load(abs_path.c_str(), &width, &height, &channels, STBI_rgb);
      if (data == nullptr) {
        throw std::runtime_error("Failed to load image: " + std::string(stbi_failure_reason()));
      }
      channels = STBI_rgb;
      img = std::vector<uint8_t>(data, data + (static_cast<ptrdiff_t>(width * height * channels)));
      stbi_image_free(data);
      if (std::cmp_not_equal(width, height)) {
        throw std::runtime_error("width != height: ");
      }
    }

    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = width - height + std::min(std::accumulate(img.begin(), img.end(), 0), channels);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (input_data_ == output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = 0;
};

namespace {

TEST_P(ZaharovGRunFuncTestsThreads, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(3, "3"), std::make_tuple(5, "5"), std::make_tuple(7, "7")};

const auto kLinContrStrsList =
    std::tuple_cat(ppc::util::AddFuncTask<ZaharovGLinContrStrALL, InType>(kTestParam, PPC_SETTINGS_zaharov_g_linear_contrast_stretch),
                   ppc::util::AddFuncTask<ZaharovGLinContrStrOMP, InType>(kTestParam, PPC_SETTINGS_zaharov_g_linear_contrast_stretch),
                   ppc::util::AddFuncTask<ZaharovGLinContrStrSEQ, InType>(kTestParam, PPC_SETTINGS_zaharov_g_linear_contrast_stretch),
                   ppc::util::AddFuncTask<ZaharovGLinContrStrSTL, InType>(kTestParam, PPC_SETTINGS_zaharov_g_linear_contrast_stretch),
                   ppc::util::AddFuncTask<ZaharovGLinContrStrTBB, InType>(kTestParam, PPC_SETTINGS_zaharov_g_linear_contrast_stretch));

const auto kGtestValues = ppc::util::ExpandToValues(kLinContrStrsList);

const auto kPerfTestName = ZaharovGRunFuncTestsThreads::PrintFuncTestName<ZaharovGRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, ZaharovGRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace zaharov_g_linear_contrast_stretch
