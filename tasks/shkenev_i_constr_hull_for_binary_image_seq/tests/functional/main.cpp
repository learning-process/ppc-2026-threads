#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "shkenev_i_constr_hull_for_binary_image_seq/common/include/common.hpp"
#include "shkenev_i_constr_hull_for_binary_image_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace shkenev_i_constr_hull_for_binary_image_seq {

namespace {

struct TestCase {
  BinaryImage image;
  std::vector<std::vector<Point>> expected;
};

BinaryImage MakeImage(int width, int height, uint8_t value = 0) {
  BinaryImage img;
  img.width = width;
  img.height = height;
  img.pixels.assign(static_cast<size_t>(width) * static_cast<size_t>(height), value);
  return img;
}

void SetPixel(BinaryImage &img, int col, int row, uint8_t value) {
  size_t idx = (static_cast<size_t>(row) * static_cast<size_t>(img.width)) + static_cast<size_t>(col);
  img.pixels[idx] = value;
}

TestCase Case1() {
  TestCase tc;
  tc.image = MakeImage(5, 5);
  SetPixel(tc.image, 2, 2, 200);
  tc.expected = {{{2, 2}}};
  return tc;
}

TestCase Case2() {
  TestCase tc;
  tc.image = MakeImage(6, 6);
  SetPixel(tc.image, 1, 1, 255);
  SetPixel(tc.image, 4, 4, 255);
  tc.expected = {{{1, 1}}, {{4, 4}}};
  return tc;
}

TestCase Case3() {
  TestCase tc;
  tc.image = MakeImage(7, 3);
  SetPixel(tc.image, 2, 1, 255);
  SetPixel(tc.image, 3, 1, 255);
  SetPixel(tc.image, 4, 1, 255);
  tc.expected = {{{2, 1}, {4, 1}}};
  return tc;
}

TestCase Case4() {
  TestCase tc;
  tc.image = MakeImage(8, 8);
  for (int row = 2; row <= 5; ++row) {
    for (int col = 3; col <= 6; ++col) {
      SetPixel(tc.image, col, row, 255);
    }
  }
  tc.expected = {{{3, 2}, {6, 2}, {6, 5}, {3, 5}}};
  return tc;
}

TestCase Case5() {
  TestCase tc;
  tc.image = MakeImage(9, 9);
  for (int row = 0; row < 9; ++row) {
    for (int col = 0; col < 9; ++col) {
      if (std::abs(col - 4) + std::abs(row - 4) <= 4) {
        SetPixel(tc.image, col, row, 255);
      }
    }
  }
  tc.expected = {{{0, 4}, {4, 0}, {8, 4}, {4, 8}}};
  return tc;
}

const std::vector<TestCase> &GetCases() {
  static std::vector<TestCase> cases = {Case1(), Case2(), Case3(), Case4(), Case5()};
  return cases;
}

const TestCase &GetCase(int id) {
  return GetCases()[static_cast<size_t>(id)];
}

std::vector<Point> Normalize(const std::vector<Point> &hull) {
  std::vector<Point> result = hull;
  std::ranges::sort(result, [](const Point &a, const Point &b) { return (a.y != b.y) ? (a.y < b.y) : (a.x < b.x); });
  auto [first, last] = std::ranges::unique(result);
  result.erase(first, last);
  return result;
}

std::vector<std::vector<Point>> NormalizeAll(const std::vector<std::vector<Point>> &hulls) {
  std::vector<std::vector<Point>> result;
  result.reserve(hulls.size());
  for (const auto &h : hulls) {
    result.push_back(Normalize(h));
  }
  std::ranges::sort(result);
  return result;
}

bool Compare(const std::vector<std::vector<Point>> &a, const std::vector<std::vector<Point>> &b) {
  return NormalizeAll(a) == NormalizeAll(b);
}

}  // namespace

class ShkenevIConstrHullFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &p) {
    return std::to_string(std::get<0>(p)) + "_" + std::get<1>(p);
  }

 protected:
  bool CheckTestOutputData(OutType &out) override {
    auto param = std::get<2>(GetParam());
    int id = std::get<0>(param);
    const auto &tc = GetCase(id);
    return Compare(tc.expected, out.convex_hulls);
  }

  InType GetTestInputData() override {
    auto param = std::get<2>(GetParam());
    int id = std::get<0>(param);
    return GetCase(id).image;
  }
};

namespace {

TEST_P(ShkenevIConstrHullFuncTests, Test) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kParams = {std::make_tuple(0, "one"), std::make_tuple(1, "two"),
                                         std::make_tuple(2, "three"), std::make_tuple(3, "four"),
                                         std::make_tuple(4, "fivee")};

const auto kTasks = ppc::util::AddFuncTask<ShkenevIConstrHullSeq, InType>(
    kParams, PPC_SETTINGS_shkenev_i_constr_hull_for_binary_image_seq);

const auto kValues = ppc::util::ExpandToValues(kTasks);
const auto kName = ShkenevIConstrHullFuncTests::PrintFuncTestName<ShkenevIConstrHullFuncTests>;

INSTANTIATE_TEST_SUITE_P(ShkenevIConstrHull, ShkenevIConstrHullFuncTests, kValues, kName);

}  // namespace

}  // namespace shkenev_i_constr_hull_for_binary_image_seq
