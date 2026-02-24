#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "romanova_v_linear_histogram_stretch/common/include/common.hpp"
#include "romanova_v_linear_histogram_stretch/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace romanova_v_linear_histogram_stretch_threads {

class RomanovaVRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    if (std::get<0>(test_param)) {
      return "loaded_pic_" + std::get<2>(test_param);
    }
    return "generated_with_size_" + std::to_string(std::get<1>(test_param));
  }

 protected:
  static std::vector<uint8_t> LoadImg(const std::string &name) {
    int width = -1;
    int height = -1;
    int channels = -1;
    std::vector<uint8_t> img;
    // Read image in RGB to ensure consistent channel count
    {
      std::string abs_path =
          ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_romanova_v_linear_histogram_stretch), name + ".png");
      auto *data = stbi_load(abs_path.c_str(), &width, &height, &channels, STBI_grey);
      if (data == nullptr) {
        throw std::runtime_error("Failed to load image: " + std::string(stbi_failure_reason()));
      }
      channels = STBI_rgb;
      img = std::vector<uint8_t>(data, data + (static_cast<ptrdiff_t>(width * height)));
      stbi_image_free(data);
    }
    return img;
  }

  static std::vector<uint8_t> MakeImg(size_t width, size_t height, uint8_t low = 75, uint8_t range = 50) {
    std::vector<uint8_t> img(width * height);
    for (size_t i = 0; i < width; i++) {
      for (size_t j = 0; j < height; j++) {
        img[(i * height) + j] = low + ((j * 13 + i * 109) % range);
      }
    }
    return img;
  }

  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    bool type = std::get<0>(params);
    if (type) {
      input_data_ = LoadImg(std::get<2>(params));
    } else {
      input_data_ = MakeImg(std::get<1>(params), std::get<1>(params));
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    auto [min_it, max_it] = std::ranges::minmax_element(output_data);
    return (*min_it == 0 && *max_it == 255) || (*min_it == *max_it);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(RomanovaVRunFuncTestsThreads, PicTests) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 4> kTestParam = {std::make_tuple(0, 3, ""), std::make_tuple(0, 5, ""),
                                            std::make_tuple(0, 70, ""), std::make_tuple(1, 0, "grey")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<RomanovaVLinHistogramStretchSEQ, InType>(
    kTestParam, PPC_SETTINGS_romanova_v_linear_histogram_stretch));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = RomanovaVRunFuncTestsThreads::PrintFuncTestName<RomanovaVRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(PicTests, RomanovaVRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace romanova_v_linear_histogram_stretch_threads
