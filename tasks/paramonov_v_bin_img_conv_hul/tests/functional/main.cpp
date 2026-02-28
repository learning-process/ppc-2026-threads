#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "paramonov_v_bin_img_conv_hul/common/include/common.hpp"
#include "paramonov_v_bin_img_conv_hul/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace paramonov_v_bin_img_conv_hul {

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
    img.pixels[static_cast<size_t>(row) * img.cols + col] = value;
  }
}

bool PointsEqual(const PixelPoint &a, const PixelPoint &b) {
  return a.row == b.row && a.col == b.col;
}

bool HullsEqual(const std::vector<PixelPoint> &h1, const std::vector<PixelPoint> &h2) {
  if (h1.size() != h2.size()) {
    return false;
  }

  std::vector<PixelPoint> sorted1 = h1;
  std::vector<PixelPoint> sorted2 = h2;

  std::sort(sorted1.begin(), sorted1.end());
  std::sort(sorted2.begin(), sorted2.end());

  return std::equal(sorted1.begin(), sorted1.end(), sorted2.begin(), PointsEqual);
}

struct TestScenario {
  GrayImage image;
  std::vector<std::vector<PixelPoint>> expected_hulls;
  std::string description;
};

std::vector<TestScenario> GenerateTestScenarios() {
  std::vector<TestScenario> scenarios;

  // Сценарий 1: Одна точка
  {
    TestScenario ts;
    ts.image = CreateTestImage(5, 5);
    SetPixel(ts.image, 2, 2);
    ts.expected_hulls = {{{2, 2}}};
    ts.description = "single_pixel";
    scenarios.push_back(ts);
  }

  // Сценарий 2: Две отдельные точки
  {
    TestScenario ts;
    ts.image = CreateTestImage(8, 8);
    SetPixel(ts.image, 1, 1);
    SetPixel(ts.image, 6, 6);
    ts.expected_hulls = {{{1, 1}}, {{6, 6}}};
    ts.description = "two_isolated_pixels";
    scenarios.push_back(ts);
  }

  // Сценарий 3: Вертикальная линия
  {
    TestScenario ts;
    ts.image = CreateTestImage(7, 7);
    for (int r = 1; r <= 5; ++r) {
      SetPixel(ts.image, r, 3);
    }
    ts.expected_hulls = {{{1, 3}, {5, 3}}};
    ts.description = "vertical_line";
    scenarios.push_back(ts);
  }

  // Сценарий 4: Горизонтальная линия
  {
    TestScenario ts;
    ts.image = CreateTestImage(7, 7);
    for (int c = 1; c <= 5; ++c) {
      SetPixel(ts.image, 3, c);
    }
    ts.expected_hulls = {{{3, 1}, {3, 5}}};
    ts.description = "horizontal_line";
    scenarios.push_back(ts);
  }

  // Сценарий 5: Прямоугольник 4x3
  {
    TestScenario ts;
    ts.image = CreateTestImage(10, 10);
    for (int r = 2; r <= 5; ++r) {
      for (int c = 3; c <= 6; ++c) {
        SetPixel(ts.image, r, c);
      }
    }
    ts.expected_hulls = {{{2, 3}, {2, 6}, {5, 6}, {5, 3}}};
    ts.description = "rectangle_4x3";
    scenarios.push_back(ts);
  }

  // Сценарий 6: Два отдельных прямоугольника
  {
    TestScenario ts;
    ts.image = CreateTestImage(15, 15);

    for (int r = 2; r <= 4; ++r) {
      for (int c = 2; c <= 4; ++c) {
        SetPixel(ts.image, r, c);
      }
    }

    for (int r = 9; r <= 11; ++r) {
      for (int c = 9; c <= 11; ++c) {
        SetPixel(ts.image, r, c);
      }
    }

    ts.expected_hulls = {{{2, 2}, {2, 4}, {4, 4}, {4, 2}}, {{9, 9}, {9, 11}, {11, 11}, {11, 9}}};
    ts.description = "two_separate_rectangles";
    scenarios.push_back(ts);
  }

  // Сценарий 7: Пустое изображение
  {
    TestScenario ts;
    ts.image = CreateTestImage(10, 10);
    ts.expected_hulls = {};
    ts.description = "empty_image";
    scenarios.push_back(ts);
  }

  return scenarios;
}

}  // namespace

class ConvexHullFunctionalTest : public ppc::util::BaseRunFuncTests<InputType, OutputType, TestScenario> {
 protected:
  bool CheckTestOutputData(OutputType &output) override {
    TestScenario param = std::get<2>(GetParam());

    if (output.size() != param.expected_hulls.size()) {
      return false;
    }

    auto hull_comparator = [](const std::vector<PixelPoint> &a, const std::vector<PixelPoint> &b) {
      if (a.empty() || b.empty()) {
        return a.size() < b.size();
      }
      return (a[0].row == b[0].row) ? a[0].col < b[0].col : a[0].row < b[0].row;
    };

    std::vector<std::vector<PixelPoint>> sorted_output = output;
    std::vector<std::vector<PixelPoint>> sorted_expected = param.expected_hulls;

    std::sort(sorted_output.begin(), sorted_output.end(), hull_comparator);
    std::sort(sorted_expected.begin(), sorted_expected.end(), hull_comparator);

    for (size_t i = 0; i < sorted_output.size(); ++i) {
      if (!HullsEqual(sorted_output[i], sorted_expected[i])) {
        return false;
      }
    }

    return true;
  }

  InputType GetTestInputData() override {
    return std::get<2>(GetParam()).image;
  }
};

TEST_P(ConvexHullFunctionalTest, Run) {
  ExecuteTest(GetParam());
}

namespace {

const auto kTestScenarios = GenerateTestScenarios();

const auto kTestTasks =
    ppc::util::AddFuncTask<ConvexHullSequential, InputType>(kTestScenarios, PPC_SETTINGS_paramonov_v_bin_img_conv_hul);

const auto kTestValues = ppc::util::ExpandToValues(kTestTasks);

INSTANTIATE_TEST_SUITE_P(ParamonovHullTests, ConvexHullFunctionalTest, kTestValues,
                         [](const testing::TestParamInfo<ConvexHullFunctionalTest::ParamType> &info) {
                           return std::get<2>(info.param).description;
                         });

}  // namespace

}  // namespace paramonov_v_bin_img_conv_hul
