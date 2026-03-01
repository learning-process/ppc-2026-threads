#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "fedoseev_linear_image_filtering_vertical/common/include/common.hpp"
#include "fedoseev_linear_image_filtering_vertical/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace fedoseev_linear_image_filtering_vertical {
namespace {

Image ReferenceFilter(const Image &input) {
  int w = input.width;
  int h = input.height;
  const std::vector<int> &src = input.data;
  std::vector<int> dst(static_cast<size_t>(w) * static_cast<size_t>(h), 0);

  const std::array<std::array<int, 3>, 3> kernel = {{{{1, 2, 1}}, {{2, 4, 2}}, {{1, 2, 1}}}};
  const int kernel_sum = 16;

  auto get = [&](int col, int row) -> int {
    col = std::clamp(col, 0, w - 1);
    row = std::clamp(row, 0, h - 1);
    return src[(static_cast<size_t>(row) * static_cast<size_t>(w)) + static_cast<size_t>(col)];
  };

  for (int row = 0; row < h; ++row) {
    for (int col = 0; col < w; ++col) {
      int sum = 0;
      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          sum += get(col + kx, row + ky) * kernel.at(ky + 1).at(kx + 1);
        }
      }
      dst[(static_cast<size_t>(row) * static_cast<size_t>(w)) + static_cast<size_t>(col)] = sum / kernel_sum;
    }
  }
  return {w, h, dst};
}

void FillConst(Image &img, [[maybe_unused]] int size) {
  std::ranges::fill(img.data, 128);
}

void FillGrad(Image &img, [[maybe_unused]] int size) {
  for (size_t i = 0; i < img.data.size(); ++i) {
    img.data[i] = static_cast<int>(i) % 256;
  }
}

void FillRand(Image &img, [[maybe_unused]] int size) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dist(0, 255);
  for (auto &v : img.data) {
    v = dist(gen);
  }
}

void FillCheckerboard(Image &img, int size) {
  const int cell = 16;
  for (int row = 0; row < size; ++row) {
    for (int col = 0; col < size; ++col) {
      size_t idx = (static_cast<size_t>(row) * static_cast<size_t>(size)) + static_cast<size_t>(col);
      img.data[idx] = (((col / cell) + (row / cell)) % 2 != 0) ? 255 : 0;
    }
  }
}

Image GenerateImage(int size, const std::string &type) {
  Image img;
  img.width = size;
  img.height = size;
  img.data.resize(static_cast<size_t>(size) * static_cast<size_t>(size));

  static const std::unordered_map<std::string, std::function<void(Image &, int)>> kFillers = {
      {"const", FillConst}, {"grad", FillGrad}, {"rand", FillRand}, {"check", FillCheckerboard}};

  auto it = kFillers.find(type);
  if (it != kFillers.end()) {
    it->second(img, size);
  } else {
    throw std::invalid_argument("Unknown type");
  }
  return img;
}

}  // namespace

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

std::array<TestType, kNumParams> GenerateParams() {
  std::array<TestType, kNumParams> params;
  auto *it = params.data();
  for (int s : kSizes) {
    for (const char *t : kTypes) {
      *it++ = std::make_tuple(s, std::string(t));
    }
  }
  return params;
}

const auto kTestParams = GenerateParams();

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<LinearImageFilteringVerticalSeq, Image>(
    kTestParams, PPC_SETTINGS_fedoseev_linear_image_filtering_vertical));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kTestName = FedoseevFuncTest::PrintFuncTestName<FedoseevFuncTest>;

INSTANTIATE_TEST_SUITE_P(ImageFilteringFuncTests, FedoseevFuncTest, kGtestValues, kTestName);

}  // namespace

}  // namespace fedoseev_linear_image_filtering_vertical
