#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "dorogin_v_bin_img_conv_hull/common/include/common.hpp"
#include "dorogin_v_bin_img_conv_hull/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace dorogin_v_bin_img_conv_hull {

namespace {

struct TestCase {
  BinaryImage image;
  std::vector<std::vector<Point>> expected;
};

BinaryImage MakeImage(int width, int height) {
  BinaryImage img;
  img.width = width;
  img.height = height;
  img.pixels.assign(static_cast<size_t>(width) * height, 0);
  return img;
}

void Set(BinaryImage &img, int x, int y) {
  img.pixels[(static_cast<size_t>(y) * img.width) + x] = 255;
}

TestCase Case0() {
  TestCase tc;
  tc.image = MakeImage(4, 4);
  Set(tc.image, 1, 2);
  tc.expected = {{{1, 2}}};
  return tc;
}

TestCase Case1() {
  TestCase tc;
  tc.image = MakeImage(7, 5);
  Set(tc.image, 0, 0);
  Set(tc.image, 6, 4);
  tc.expected = {{{0, 0}}, {{6, 4}}};
  return tc;
}

TestCase Case2() {
  TestCase tc;
  tc.image = MakeImage(5, 6);
  for (int col = 1; col <= 4; ++col) {
    Set(tc.image, 3, col);
  }
  tc.expected = {{{3, 1}, {3, 4}}};
  return tc;
}

TestCase Case3() {
  TestCase tc;
  tc.image = MakeImage(10, 6);
  for (int col = 1; col <= 4; ++col) {
    for (int row = 2; row <= 7; ++row) {
      Set(tc.image, row, col);
    }
  }

  tc.expected = {{{2, 1}, {7, 1}, {7, 4}, {2, 4}}};
  return tc;
}

TestCase Case4() {
  TestCase tc;
  tc.image = MakeImage(11, 11);

  for (int col = 0; col < 11; ++col) {
    for (int row = 0; row < 11; ++row) {
      if (std::abs(row - 5) + std::abs(col - 5) <= 5) {
        Set(tc.image, row, col);
      }
    }
  }

  tc.expected = {{{0, 5}, {5, 0}, {10, 5}, {5, 10}}};
  return tc;
}

const std::vector<TestCase> &GetCases() {
  static std::vector<TestCase> cases = {Case0(), Case1(), Case2(), Case3(), Case4()};
  return cases;
}

const TestCase &GetCase(int id) {
  return GetCases()[static_cast<size_t>(id)];
}

std::vector<Point> Normalize(const std::vector<Point> &hull) {
  std::vector<Point> result = hull;
  std::sort(result.begin(), result.end(),
            [](const Point &a, const Point &b) { return (a.y == b.y) ? a.x < b.x : a.y < b.y; });
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
}

std::vector<std::vector<Point>> NormalizeAll(const std::vector<std::vector<Point>> &hulls) {
  std::vector<std::vector<Point>> result;
  result.reserve(hulls.size());
  for (const auto &h : hulls) {
    result.push_back(Normalize(h));
  }

  std::sort(result.begin(), result.end());
  return result;
}

bool Compare(const std::vector<std::vector<Point>> &a, const std::vector<std::vector<Point>> &b) {
  return NormalizeAll(a) == NormalizeAll(b);
}

}  // namespace

class DorogiVBinImgConvHullFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &p) {
    return std::to_string(std::get<0>(p)) + "_" + std::get<1>(p);
  }

 protected:
  bool CheckTestOutputData(OutType &out) override {
    auto param = std::get<2>(GetParam());
    int id = std::get<0>(param);
    return Compare(GetCase(id).expected, out.convex_hulls);
  }

  InType GetTestInputData() override {
    auto param = std::get<2>(GetParam());
    return GetCase(std::get<0>(param)).image;
  }
};

TEST_P(DorogiVBinImgConvHullFuncTests, Run) {
  ExecuteTest(GetParam());
}

namespace {

const std::array<TestType, 5> kParams = {std::make_tuple(0, "single"), std::make_tuple(1, "two_points"),
                                         std::make_tuple(2, "vertical_line"), std::make_tuple(3, "rectangle"),
                                         std::make_tuple(4, "diamond_large")};

const auto kTasks =
    ppc::util::AddFuncTask<DoroginVBinImgConvHullSeq, InType>(kParams, PPC_SETTINGS_dorogin_v_bin_img_conv_hull);

const auto kValues = ppc::util::ExpandToValues(kTasks);

INSTANTIATE_TEST_SUITE_P(DHull, DorogiVBinImgConvHullFuncTests, kValues,
                         DorogiVBinImgConvHullFuncTests::PrintFuncTestName<DorogiVBinImgConvHullFuncTests>);

}  // namespace

}  // namespace dorogin_v_bin_img_conv_hull
