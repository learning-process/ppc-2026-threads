#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "egorova_l_binary_convex_hull/common/include/common.hpp"
#include "egorova_l_binary_convex_hull/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace egorova_l_binary_convex_hull {

using TestType = std::tuple<InType, std::vector<std::vector<Point>>, std::string>;

// Helper functions moved outside the class to reduce complexity
namespace {
bool ArePointsEqual(const Point &p1, const Point &p2) {
  return p1.x == p2.x && p1.y == p2.y;
}

void SortPoints(std::vector<Point> &points) {
  std::ranges::sort(points,
                    [](const Point &lhs, const Point &rhs) { return std::tie(lhs.x, lhs.y) < std::tie(rhs.x, rhs.y); });
}

bool AreHullsEqual(const std::vector<Point> &hull1, const std::vector<Point> &hull2) {
  if (hull1.size() != hull2.size()) {
    return false;
  }

  std::vector<Point> sorted1 = hull1;
  std::vector<Point> sorted2 = hull2;

  SortPoints(sorted1);
  SortPoints(sorted2);

  for (size_t i = 0; i < sorted1.size(); ++i) {
    if (!ArePointsEqual(sorted1[i], sorted2[i])) {
      return false;
    }
  }
  return true;
}

void SortHulls(std::vector<std::vector<Point>> &hulls) {
  for (auto &hull : hulls) {
    SortPoints(hull);
  }

  std::ranges::sort(hulls, [](const std::vector<Point> &a, const std::vector<Point> &b) {
    if (a.empty() || b.empty()) {
      return a.size() < b.size();
    }
    return std::tie(a[0].x, a[0].y) < std::tie(b[0].x, b[0].y);
  });
}

}  // namespace

class EgorovaLFuncTest : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_result_ = std::get<1>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_result_.size()) {
      return false;
    }

    std::vector<std::vector<Point>> sorted_output = output_data;
    std::vector<std::vector<Point>> sorted_expected = expected_result_;

    SortHulls(sorted_output);
    SortHulls(sorted_expected);

    for (size_t i = 0; i < sorted_output.size(); ++i) {
      if (!AreHullsEqual(sorted_output[i], sorted_expected[i])) {
        return false;
      }
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_result_;
};

namespace {

InType CreateEmptyImage(int width, int height) {
  InType img;
  img.width = width;
  img.height = height;
  img.data.assign(static_cast<size_t>(width) * static_cast<size_t>(height), 0);
  return img;
}

void DrawRectangle(InType &img, int x1, int y1, int x2, int y2) {
  for (int row = y1; row <= y2; ++row) {
    for (int col = x1; col <= x2; ++col) {
      if (col >= 0 && col < img.width && row >= 0 && row < img.height) {
        img.data[(static_cast<size_t>(row) * static_cast<size_t>(img.width)) + static_cast<size_t>(col)] = 255;
      }
    }
  }
}

std::vector<Point> GetRectangleHull(int x1, int y1, int x2, int y2) {
  return {{x1, y1}, {x2, y1}, {x2, y2}, {x1, y2}};
}

const std::array<TestType, 7> kTestParams = {{
    // Test 1: Single 3x3 square
    std::make_tuple(
        []() {
          auto img = CreateEmptyImage(10, 10);
          DrawRectangle(img, 1, 1, 3, 3);
          return img;
        }(),
        std::vector<std::vector<Point>>{GetRectangleHull(1, 1, 3, 3)},
        "single_square"),

    // Test 2: Two separate squares
    std::make_tuple(
        []() {
          auto img = CreateEmptyImage(20, 20);
          DrawRectangle(img, 2, 2, 4, 4);
          DrawRectangle(img, 10, 10, 12, 12);
          return img;
        }(),
        std::vector<std::vector<Point>>{GetRectangleHull(2, 2, 4, 4),
                                        GetRectangleHull(10, 10, 12, 12)},
        "two_squares"),

    // Test 3: 5x3 rectangle
    std::make_tuple(
        []() {
          auto img = CreateEmptyImage(15, 15);
          DrawRectangle(img, 2, 3, 6, 5);
          return img;
        }(),
        std::vector<std::vector<Point>>{GetRectangleHull(2, 3, 6, 5)},
        "rectangle"),

    // Test 4: Three components
    std::make_tuple(
        []() {
          auto img = CreateEmptyImage(30, 30);
          DrawRectangle(img, 1, 1, 3, 3);
          DrawRectangle(img, 10, 10, 12, 12);
          DrawRectangle(img, 20, 5, 22, 7);
          return img;
        }(),
        std::vector<std::vector<Point>>{GetRectangleHull(1, 1, 3, 3),
                                        GetRectangleHull(10, 10, 12, 12),
                                        GetRectangleHull(20, 5, 22, 7)},
        "three_components"),

    // Test 5: Empty image
    std::make_tuple(CreateEmptyImage(10, 10), std::vector<std::vector<Point>>{},
                    "empty_image"),

    // Test 6: Horizontal line (should give 2 points)
    std::make_tuple(
        []() {
          auto img = CreateEmptyImage(10, 10);
          for (int col = 2; col <= 5; ++col) {
            img.data[(static_cast<size_t>(5) * static_cast<size_t>(img.width)) + 
                      static_cast<size_t>(col)] = 255;
          }
          return img;
        }(),
        std::vector<std::vector<Point>>{{{2, 5}, {5, 5}}},
        "horizontal_line"),

    // Test 7: Vertical line (should give 2 points)
    std::make_tuple(
        []() {
          auto img = CreateEmptyImage(10, 10);
          for (int row = 2; row <= 5; ++row) {
            img.data[(static_cast<size_t>(row) * static_cast<size_t>(img.width)) + 
                      static_cast<size_t>(5)] = 255;
          }
          return img;
        }(),
        std::vector<std::vector<Point>>{{{5, 2}, {5, 5}}},
        "vertical_line")}};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<BinaryConvexHullSEQ, InType>(kTestParams, PPC_SETTINGS_egorova_l_binary_convex_hull));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kFuncTestName = EgorovaLFuncTest::PrintFuncTestName<EgorovaLFuncTest>;

INSTANTIATE_TEST_SUITE_P(BinaryConvexHullTests, EgorovaLFuncTest, kGtestValues, kFuncTestName);

TEST_P(EgorovaLFuncTest, RunFunctionalTests) {
  ExecuteTest(GetParam());
}

}  // namespace

}  // namespace egorova_l_binary_convex_hull
