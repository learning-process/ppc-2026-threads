#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "kamalagin_a_binary_image_convex_hull/common/include/common.hpp"
#include "kamalagin_a_binary_image_convex_hull/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kamalagin_a_binary_image_convex_hull {

namespace {

bool PointLess(const Point &a, const Point &b) {
  return a.x != b.x ? a.x < b.x : a.y < b.y;
}

void NormalizeHull(Hull &h) {
  std::ranges::sort(h, PointLess);
}

void NormalizeHullList(HullList &list) {
  for (auto &h : list) {
    NormalizeHull(h);
  }
  std::ranges::sort(list, [](const Hull &a, const Hull &b) {
    if (a.size() != b.size()) {
      return a.size() < b.size();
    }
    for (size_t i = 0; i < a.size(); ++i) {
      if (a[i].x != b[i].x) {
        return a[i].x < b[i].x;
      }
      if (a[i].y != b[i].y) {
        return a[i].y < b[i].y;
      }
    }
    return false;
  });
}

BinaryImage MakeEmptyImage() {
  BinaryImage img;
  img.rows = 0;
  img.cols = 0;
  return img;
}

BinaryImage MakeSinglePixel() {
  BinaryImage img;
  img.rows = 1;
  img.cols = 1;
  img.data = {1};
  return img;
}

BinaryImage MakeTwoIsolatedPixels() {
  BinaryImage img;
  img.rows = 3;
  img.cols = 3;
  img.data.assign(9, 0);
  img.data[0] = 1;
  img.data[8] = 1;
  return img;
}

BinaryImage MakeVerticalLine() {
  BinaryImage img;
  img.rows = 5;
  img.cols = 1;
  img.data.assign(5, 1);
  return img;
}

BinaryImage MakeHorizontalLine() {
  BinaryImage img;
  img.rows = 1;
  img.cols = 5;
  img.data.assign(5, 1);
  return img;
}

BinaryImage MakeFilledRectangle() {
  BinaryImage img;
  img.rows = 2;
  img.cols = 3;
  img.data.assign(6, 1);
  return img;
}

BinaryImage MakeTwoComponents() {
  BinaryImage img;
  img.rows = 1;
  img.cols = 5;
  img.data = {1, 0, 0, 0, 1};
  return img;
}

HullList ExpectedEmpty() {
  return {};
}

HullList ExpectedSinglePixel() {
  return {{Point{.x = 0, .y = 0}}};
}

HullList ExpectedTwoIsolatedPixels() {
  return {{Point{.x = 0, .y = 0}}, {Point{.x = 2, .y = 2}}};
}

HullList ExpectedVerticalLine() {
  return {{Point{.x = 0, .y = 0}, Point{.x = 0, .y = 4}}};
}

HullList ExpectedHorizontalLine() {
  return {{Point{.x = 0, .y = 0}, Point{.x = 4, .y = 0}}};
}

HullList ExpectedFilledRectangle() {
  return {{Point{.x = 0, .y = 0}, Point{.x = 0, .y = 1}, Point{.x = 2, .y = 0}, Point{.x = 2, .y = 1}}};
}

HullList ExpectedTwoComponents() {
  return {{Point{.x = 0, .y = 0}}, {Point{.x = 4, .y = 0}}};
}

using MakeFn = BinaryImage (*)();
using ExpectedFn = HullList (*)();

struct TestCaseEntry {
  const char *name;
  MakeFn make;
  ExpectedFn expected;
};

const std::array<TestCaseEntry, 8> kTestCaseTable = {{
    {.name = "empty", .make = MakeEmptyImage, .expected = ExpectedEmpty},
    {.name = "single_pixel", .make = MakeSinglePixel, .expected = ExpectedSinglePixel},
    {.name = "two_isolated", .make = MakeTwoIsolatedPixels, .expected = ExpectedTwoIsolatedPixels},
    {.name = "vertical_line", .make = MakeVerticalLine, .expected = ExpectedVerticalLine},
    {.name = "horizontal_line", .make = MakeHorizontalLine, .expected = ExpectedHorizontalLine},
    {.name = "filled_rect", .make = MakeFilledRectangle, .expected = ExpectedFilledRectangle},
    {.name = "two_components", .make = MakeTwoComponents, .expected = ExpectedTwoComponents},
}};

std::pair<BinaryImage, HullList> GetTestDataForName(const std::string &name) {
  for (const auto &entry : kTestCaseTable) {
    if (name == entry.name) {
      return {entry.make(), entry.expected()};
    }
  }
  return {MakeEmptyImage(), ExpectedEmpty()};
}

}  // namespace

class KamalaginRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    const std::string &name = std::get<1>(params);
    std::tie(input_data_, expected_) = GetTestDataForName(name);
    NormalizeHullList(expected_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    NormalizeHullList(output_data);
    if (output_data.size() != expected_.size()) {
      return false;
    }
    for (size_t i = 0; i < expected_.size(); ++i) {
      if (output_data[i].size() != expected_[i].size()) {
        return false;
      }
      for (size_t j = 0; j < expected_[i].size(); ++j) {
        if (output_data[i][j].x != expected_[i][j].x || output_data[i][j].y != expected_[i][j].y) {
          return false;
        }
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  HullList expected_;
};

namespace {

TEST_P(KamalaginRunFuncTests, BinaryImageConvexHull) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 7> kTestParam = {
    std::make_tuple(0, "empty"),          std::make_tuple(1, "single_pixel"),    std::make_tuple(2, "two_isolated"),
    std::make_tuple(3, "vertical_line"),  std::make_tuple(4, "horizontal_line"), std::make_tuple(5, "filled_rect"),
    std::make_tuple(6, "two_components"),
};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<KamalaginABinaryImageConvexHullSEQ, InType>(
    kTestParam, PPC_SETTINGS_kamalagin_a_binary_image_convex_hull));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kPerfTestName = KamalaginRunFuncTests::PrintFuncTestName<KamalaginRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(ConvexHullTests, KamalaginRunFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kamalagin_a_binary_image_convex_hull
