#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "paramonov_v_bin_img_conv_hul/common/include/common.hpp"
#include "paramonov_v_bin_img_conv_hul/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
// GTestParamIndex используется через func_test_util.hpp, поэтому отдельный include не требуется

namespace paramonov_v_bin_img_conv_hul {

using TestCase = std::tuple<GrayImage, std::vector<std::vector<PixelPoint>>, std::string>;

namespace {

GrayImage CreateTestImage(int rows, int cols) {
  GrayImage img;
  img.rows = rows;
  img.cols = cols;
  img.pixels.assign(static_cast<size_t>(rows) * cols, 0);
  return img;
}

void SetPixel(GrayImage &img, int row, int col, uint8_t value = 255) {
  if (row >= 0 && row < img.rows && col >= 0 && col < img.cols) {
    img.pixels[(static_cast<size_t>(row) * img.cols) + col] = value;
  }
}

void DrawRectangle(GrayImage &img, int start_row, int start_col, int end_row, int end_col) {
  for (int row = start_row; row <= end_row; ++row) {
    for (int col = start_col; col <= end_col; ++col) {
      SetPixel(img, row, col);
    }
  }
}

std::vector<PixelPoint> GetRectangleHull(int start_row, int start_col, int end_row, int end_col) {
  return {{start_row, start_col}, {start_row, end_col}, {end_row, end_col}, {end_row, start_col}};
}

bool PointsEqual(const PixelPoint &a, const PixelPoint &b) {
  return a.row == b.row && a.col == b.col;
}

bool ComparePoints(const PixelPoint &a, const PixelPoint &b) {
  if (a.row != b.row) {
    return a.row < b.row;
  }
  return a.col < b.col;
}

void SortPoints(std::vector<PixelPoint> &points) {
  std::ranges::sort(points, ComparePoints);
}

bool HullsEqual(const std::vector<PixelPoint> &hull1, const std::vector<PixelPoint> &hull2) {
  if (hull1.size() != hull2.size()) {
    return false;
  }

  std::vector<PixelPoint> sorted1 = hull1;
  std::vector<PixelPoint> sorted2 = hull2;

  SortPoints(sorted1);
  SortPoints(sorted2);

  for (size_t i = 0; i < sorted1.size(); ++i) {
    if (!PointsEqual(sorted1[i], sorted2[i])) {
      return false;
    }
  }
  return true;
}

bool CompareHulls(const std::vector<PixelPoint> &a, const std::vector<PixelPoint> &b) {
  if (a.empty() || b.empty()) {
    return a.size() < b.size();
  }
  if (a[0].row != b[0].row) {
    return a[0].row < b[0].row;
  }
  return a[0].col < b[0].col;
}

void SortHulls(std::vector<std::vector<PixelPoint>> &hulls) {
  for (auto &hull : hulls) {
    SortPoints(hull);
  }

  std::ranges::sort(hulls, CompareHulls);
}

}  // namespace

class ConvexHullFuncTest : public ppc::util::BaseRunFuncTests<InputType, OutputType, TestCase> {
 public:
  static std::string PrintTestParam(const TestCase &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    // Используем прямой доступ к параметрам через индексы, чтобы избежать предупреждения
    const auto &param_tuple = GetParam();
    const auto &test_params = std::get<2>(param_tuple);
    input_image_ = std::get<0>(test_params);
    expected_hulls_ = std::get<1>(test_params);
  }

  bool CheckTestOutputData(OutputType &output) override {
    if (output.size() != expected_hulls_.size()) {
      return false;
    }

    std::vector<std::vector<PixelPoint>> sorted_output = output;
    std::vector<std::vector<PixelPoint>> sorted_expected = expected_hulls_;

    SortHulls(sorted_output);
    SortHulls(sorted_expected);

    for (size_t i = 0; i < sorted_output.size(); ++i) {
      if (!HullsEqual(sorted_output[i], sorted_expected[i])) {
        return false;
      }
    }

    return true;
  }

  InputType GetTestInputData() override {
    return input_image_;
  }

 private:
  InputType input_image_;
  OutputType expected_hulls_;
};

namespace {

const std::array<TestCase, 8> kTestCases = {{
    // Тест 1: Одна точка
    std::make_tuple(
        []() {
          auto img = CreateTestImage(5, 5);
          SetPixel(img, 2, 2);
          return img;
        }(),
        std::vector<std::vector<PixelPoint>>{{{2, 2}}},
        "single_pixel"),

    // Тест 2: Две отдельные точки
    std::make_tuple(
        []() {
          auto img = CreateTestImage(8, 8);
          SetPixel(img, 1, 1);
          SetPixel(img, 6, 6);
          return img;
        }(),
        std::vector<std::vector<PixelPoint>>{{{1, 1}}, {{6, 6}}},
        "two_isolated_pixels"),

    // Тест 3: Вертикальная линия
    std::make_tuple(
        []() {
          auto img = CreateTestImage(7, 7);
          for (int row = 1; row <= 5; ++row) {
            SetPixel(img, row, 3);
          }
          return img;
        }(),
        std::vector<std::vector<PixelPoint>>{{{1, 3}, {5, 3}}},
        "vertical_line"),

    // Тест 4: Горизонтальная линия
    std::make_tuple(
        []() {
          auto img = CreateTestImage(7, 7);
          for (int col = 1; col <= 5; ++col) {
            SetPixel(img, 3, col);
          }
          return img;
        }(),
        std::vector<std::vector<PixelPoint>>{{{3, 1}, {3, 5}}},
        "horizontal_line"),

    // Тест 5: Прямоугольник
    std::make_tuple(
        []() {
          auto img = CreateTestImage(10, 10);
          DrawRectangle(img, 2, 3, 5, 6);
          return img;
        }(),
        std::vector<std::vector<PixelPoint>>{GetRectangleHull(2, 3, 5, 6)},
        "rectangle"),

    // Тест 6: Два прямоугольника
    std::make_tuple(
        []() {
          auto img = CreateTestImage(15, 15);
          DrawRectangle(img, 2, 2, 4, 4);
          DrawRectangle(img, 9, 9, 11, 11);
          return img;
        }(),
        std::vector<std::vector<PixelPoint>>{
          GetRectangleHull(2, 2, 4, 4),
          GetRectangleHull(9, 9, 11, 11)
        },
        "two_rectangles"),

    // Тест 7: Три компоненты
    std::make_tuple(
        []() {
          auto img = CreateTestImage(30, 30);
          DrawRectangle(img, 1, 1, 3, 3);
          DrawRectangle(img, 10, 10, 12, 12);
          DrawRectangle(img, 20, 5, 22, 7);
          return img;
        }(),
        std::vector<std::vector<PixelPoint>>{
          GetRectangleHull(1, 1, 3, 3),
          GetRectangleHull(10, 10, 12, 12),
          GetRectangleHull(20, 5, 22, 7)
        },
        "three_components"),

    // Тест 8: Пустое изображение
    std::make_tuple(
        CreateTestImage(10, 10),
        std::vector<std::vector<PixelPoint>>{},
        "empty_image")
}};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<ConvexHullSequential, InputType>(kTestCases, PPC_SETTINGS_paramonov_v_bin_img_conv_hul));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kFuncTestName = ConvexHullFuncTest::PrintFuncTestName<ConvexHullFuncTest>;

INSTANTIATE_TEST_SUITE_P(ParamonovHullTests, ConvexHullFuncTest, kGtestValues, kFuncTestName);

TEST_P(ConvexHullFuncTest, RunFunctionalTests) {
  ExecuteTest(GetParam());
}

}  // namespace

}  // namespace paramonov_v_bin_img_conv_hul
