#include <gtest/gtest.h>

#include <array>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "fedoseev_linear_image_filtering_vertical/common/include/common.hpp"
#include "fedoseev_linear_image_filtering_vertical/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace fedoseev_linear_image_filtering_vertical {

static Image ReferenceFilter(const Image &input) {
  int w = input.width;
  int h = input.height;
  const std::vector<int> &src = input.data;
  std::vector<int> dst(w * h, 0);

  const int kernel[3][3] = {{1, 2, 1}, {2, 4, 2}, {1, 2, 1}};
  const int kernel_sum = 16;

  auto get = [&](int x, int y) -> int {
    x = std::clamp(x, 0, w - 1);
    y = std::clamp(y, 0, h - 1);
    return src[y * w + x];
  };

  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      int sum = 0;
      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          sum += get(x + kx, y + ky) * kernel[ky + 1][kx + 1];
        }
      }
      dst[y * w + x] = sum / kernel_sum;
    }
  }
  return {w, h, dst};
}

static Image GenerateImage(int size, const std::string &type) {
  Image img;
  img.width = size;
  img.height = size;
  img.data.resize(size * size);

  if (type == "const") {
    std::fill(img.data.begin(), img.data.end(), 128);
  } else if (type == "grad") {
    for (size_t i = 0; i < img.data.size(); ++i) {
      img.data[i] = static_cast<int>(i) % 256;
    }
  } else if (type == "rand") {
    std::mt19937 gen(42);
    std::uniform_int_distribution<int> dist(0, 255);
    for (auto &v : img.data) {
      v = dist(gen);
    }
  } else if (type == "check") {
    int cell = 16;
    for (int y = 0; y < size; ++y) {
      for (int x = 0; x < size; ++x) {
        img.data[y * size + x] = ((x / cell + y / cell) % 2) ? 255 : 0;
      }
    }
  } else {
    throw std::invalid_argument("Unknown type");
  }
  return img;
}

class FedoseevFuncTest : public ppc::util::BaseRunFuncTests<Image, Image, TestType> {
 public:
  static std::string PrintTestParam(const TestType &param) {
    return std::to_string(std::get<0>(param)) + "_" + std::get<1>(param);
  }

 protected:
  void SetUp() override {
    auto param = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int size = std::get<0>(param);
    std::string type = std::get<1>(param);

    input_ = GenerateImage(size, type);
    expected_ = ReferenceFilter(input_);
  }

  bool CheckTestOutputData(Image &output_data) override {
    if (output_data.width != expected_.width || output_data.height != expected_.height) {
      return false;
    }
    return output_data.data == expected_.data;
  }

  Image GetTestInputData() override {
    return input_;
  }

 private:
  Image input_;
  Image expected_;
};

namespace {

constexpr std::array<int, 5> kSizes = {3, 5, 7, 10, 16};
constexpr std::array<const char *, 4> kTypes = {"const", "grad", "rand", "check"};

constexpr size_t kNumParams = kSizes.size() * kTypes.size();

template <size_t... Is>
constexpr std::array<TestType, kNumParams> GenerateParamsImpl(std::index_sequence<Is...>) {
  return {std::make_tuple(kSizes[Is / kTypes.size()], kTypes[Is % kTypes.size()])...};
}

constexpr auto kTestParams = GenerateParamsImpl(std::make_index_sequence<kNumParams>{});

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<LinearImageFilteringVerticalSeq, Image>(
    kTestParams, PPC_SETTINGS_fedoseev_linear_image_filtering_vertical));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kTestName = FedoseevFuncTest::PrintFuncTestName<FedoseevFuncTest>;

INSTANTIATE_TEST_SUITE_P(ImageFilteringFuncTests, FedoseevFuncTest, kGtestValues, kTestName);

}  // namespace

}  // namespace fedoseev_linear_image_filtering_vertical
