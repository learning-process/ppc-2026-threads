#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <ranges>
#include <string>
#include <tuple>
#include <vector>

#include "dorogin_v_bin_img_conv_hull_TBB/common/include/common.hpp"
#include "dorogin_v_bin_img_conv_hull_TBB/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace dorogin_v_bin_img_conv_hull_tbb {
using TestCase = std::tuple<GrayImage, std::vector<std::vector<PixelPoint>>, std::string>;

namespace {
GrayImage CreateImage(int rows, int cols) {
  GrayImage image;
  image.rows = rows;
  image.cols = cols;
  image.pixels.assign(static_cast<size_t>(rows) * static_cast<size_t>(cols), 0);
  return image;
}

void SetPixel(GrayImage &image, int row, int col, uint8_t value = 255) {
  if (row < 0 || col < 0 || row >= image.rows || col >= image.cols) {
    return;
  }
  image.pixels[(static_cast<size_t>(row) * static_cast<size_t>(image.cols)) + static_cast<size_t>(col)] = value;
}

void FillRectangle(GrayImage &image, int row0, int col0, int row1, int col1) {
  for (int row_idx = row0; row_idx <= row1; ++row_idx) {
    for (int col_idx = col0; col_idx <= col1; ++col_idx) {
      SetPixel(image, row_idx, col_idx, 255);
    }
  }
}

std::vector<PixelPoint> RectangleHull(int row0, int col0, int row1, int col1) {
  return {{row0, col0}, {row0, col1}, {row1, col1}, {row1, col0}};
}

bool PointLess(const PixelPoint &a, const PixelPoint &b) {
  if (a.row != b.row) {
    return a.row < b.row;
  }
  return a.col < b.col;
}

bool HullLess(const std::vector<PixelPoint> &a, const std::vector<PixelPoint> &b) {
  if (a.empty() || b.empty()) {
    return a.size() < b.size();
  }
  return PointLess(a.front(), b.front());
}

void SortHull(std::vector<PixelPoint> *hull) {
  std::ranges::sort(*hull, PointLess);
}

void Normalize(std::vector<std::vector<PixelPoint>> *hulls) {
  for (auto &hull : *hulls) {
    SortHull(&hull);
  }
  std::ranges::sort(*hulls, HullLess);
}
}  // namespace

class DoroginVImgConvHullTBBFuncTest : public ppc::util::BaseRunFuncTests<InputType, OutputType, TestCase> {
 public:
  static std::string PrintTestParam(const TestCase &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    const auto &input = std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_ = std::get<0>(input);
    expected_ = std::get<1>(input);
  }

  bool CheckTestOutputData(OutputType &output) override {
    if (output.size() != expected_.size()) {
      return false;
    }

    auto output_normalized = output;
    auto expected_normalized = expected_;
    Normalize(&output_normalized);
    Normalize(&expected_normalized);
    return output_normalized == expected_normalized;
  }

  InputType GetTestInputData() override {
    return input_;
  }

 private:
  InputType input_;
  OutputType expected_;
};

namespace {
const std::array<TestCase, 8> kCases = {
    std::make_tuple(
        [] {
  auto image = CreateImage(5, 5);
  SetPixel(image, 2, 2);
  return image;
}(), std::vector<std::vector<PixelPoint>>{{{2, 2}}}, "single_pixel"),
    std::make_tuple(
        [] {
  auto image = CreateImage(8, 8);
  SetPixel(image, 1, 1);
  SetPixel(image, 6, 6);
  return image;
}(), std::vector<std::vector<PixelPoint>>{{{1, 1}}, {{6, 6}}}, "two_isolated_pixels"),
    std::make_tuple(
        [] {
  auto image = CreateImage(7, 7);
  for (int row_idx = 1; row_idx <= 5; ++row_idx) {
    SetPixel(image, row_idx, 3);
  }
  return image;
}(), std::vector<std::vector<PixelPoint>>{{{1, 3}, {5, 3}}}, "vertical_line"),
    std::make_tuple(
        [] {
  auto image = CreateImage(7, 7);
  for (int col_idx = 1; col_idx <= 5; ++col_idx) {
    SetPixel(image, 3, col_idx);
  }
  return image;
}(), std::vector<std::vector<PixelPoint>>{{{3, 1}, {3, 5}}}, "horizontal_line"),
    std::make_tuple(
        [] {
  auto image = CreateImage(10, 10);
  FillRectangle(image, 2, 3, 5, 6);
  return image;
}(), std::vector<std::vector<PixelPoint>>{RectangleHull(2, 3, 5, 6)}, "rectangle"),
    std::make_tuple(
        [] {
  auto image = CreateImage(15, 15);
  FillRectangle(image, 2, 2, 4, 4);
  FillRectangle(image, 9, 9, 11, 11);
  return image;
}(), std::vector<std::vector<PixelPoint>>{RectangleHull(2, 2, 4, 4), RectangleHull(9, 9, 11, 11)}, "two_rectangles"),
    std::make_tuple(
        [] {
  auto image = CreateImage(30, 30);
  FillRectangle(image, 1, 1, 3, 3);
  FillRectangle(image, 10, 10, 12, 12);
  FillRectangle(image, 20, 5, 22, 7);
  return image;
}(),
        std::vector<std::vector<PixelPoint>>{RectangleHull(1, 1, 3, 3), RectangleHull(10, 10, 12, 12),
                                             RectangleHull(20, 5, 22, 7)},
        "three_components"),
    std::make_tuple(CreateImage(10, 10), std::vector<std::vector<PixelPoint>>{}, "empty_image")};

const auto kTasks = std::tuple_cat(
    ppc::util::AddFuncTask<DoroginVImgConvHullTBB, InputType>(kCases, PPC_SETTINGS_dorogin_v_bin_img_conv_hull_TBB));
const auto kValues = ppc::util::ExpandToValues(kTasks);
const auto kName = DoroginVImgConvHullTBBFuncTest::PrintFuncTestName<DoroginVImgConvHullTBBFuncTest>;

INSTANTIATE_TEST_SUITE_P(DoroginVBinImgConvHullTBBFunctional, DoroginVImgConvHullTBBFuncTest, kValues, kName);

TEST_P(DoroginVImgConvHullTBBFuncTest, DoroginVFuncBinaryImageConvexHullTBB) {
  ExecuteTest(GetParam());
}
}  // namespace

}  // namespace dorogin_v_bin_img_conv_hull_tbb
