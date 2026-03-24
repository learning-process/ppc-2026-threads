#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "rychkova_gauss/common/include/common.hpp"
#include "rychkova_gauss/seq/include/ops_seq.hpp"
#include "rychkova_gauss/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace rychkova_gauss {

class RychkovaGaussFunc : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return test_param;
  }

 protected:
  static Image ReadImage(const std::string &file) {
    int width = -1;
    int height = -1;
    int channels = -1;
    std::vector<uint8_t> img;
    {
      std::string abs_path =
          ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_rychkova_gauss), file);  // путь до папки дата
      auto *data = stbi_load(abs_path.c_str(), &width, &height, &channels, STBI_rgb);
      if (data == nullptr) {
        throw std::runtime_error("Failed to load image: " + std::string(stbi_failure_reason()));
      }
      channels = STBI_rgb;
      img = std::vector<uint8_t>(data,
                                 data + (static_cast<ptrdiff_t>(width * height * channels)));  // колво элем в массиве
      stbi_image_free(data);  // освобождение указателя (память)
    }
    Image inp(height, std::vector<Pixel>(width, {0, 0, 0}));
    for (int i = 0; i < width * height * channels; i++) {
      int row = (i / 3) / width;
      int column = (i / 3) % width;
      int ch = i % 3;
      if (ch == 0) {
        inp[row][column].R = img[i];
      } else if (ch == 1) {
        inp[row][column].G = img[i];
      } else if (ch == 2) {
        inp[row][column].B = img[i];
      }
    }
    return inp;
  }

  void SetUp() override {
    TestType param = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = ReadImage("input/" + param + ".png");
    expected_data_ = ReadImage("expected/" + param + "_exp.png");
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::ranges::equal(output_data, expected_data_, [](auto first, auto second) {
      return std::ranges::equal(first, second, [](auto a, auto b) {
        return (std::fabs(static_cast<int>(a.R) - static_cast<int>(b.R)) < kEPS) &&
               (std::fabs(static_cast<int>(a.G) - static_cast<int>(b.G)) < kEPS) &&
               (std::fabs(static_cast<int>(a.B) - static_cast<int>(b.B)) < kEPS);
      });
    });
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_data_;
};

namespace {

TEST_P(RychkovaGaussFunc, GaussBlur) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 8> kTestParam = {"high_freq_noize", "img_1", "img_2",     "img_3",
                                            "img_4",           "lines", "rgb_lines", "low_freq_noize"};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<RychkovaGaussSEQ, InType>(kTestParam, PPC_SETTINGS_rychkova_gauss),
                   ppc::util::AddFuncTask<RychkovaGaussTBB, InType>(kTestParam, PPC_SETTINGS_rychkova_gauss));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = RychkovaGaussFunc::PrintFuncTestName<RychkovaGaussFunc>;

INSTANTIATE_TEST_SUITE_P(RychkovaGaussFuncSuite, RychkovaGaussFunc, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace rychkova_gauss
