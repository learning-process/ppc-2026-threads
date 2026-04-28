#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>

#include "kamalagin_a_binary_image_convex_hull_omp/common/include/common.hpp"
#include "kamalagin_a_binary_image_convex_hull_omp/omp/include/ops_omp.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kamalagin_a_binary_image_convex_hull_omp {

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

struct TestCaseEntry {
  const char *name;
  BinaryImage (*make)();
  HullList expected;
};

BinaryImage MakeEmpty() {
  return BinaryImage{.rows = 0, .cols = 0, .data = {}};
}

BinaryImage MakeSinglePixel() {
  return BinaryImage{.rows = 1, .cols = 1, .data = {1}};
}

BinaryImage MakeTwoIsolated() {
  return BinaryImage{.rows = 2, .cols = 2, .data = {1, 0, 0, 1}};
}

BinaryImage MakeVerticalLine() {
  return BinaryImage{.rows = 3, .cols = 1, .data = {1, 1, 1}};
}

BinaryImage MakeHorizontalLine() {
  return BinaryImage{.rows = 1, .cols = 3, .data = {1, 1, 1}};
}

BinaryImage MakeFilledRect() {
  return BinaryImage{.rows = 2, .cols = 3, .data = {1, 1, 1, 1, 1, 1}};
}

BinaryImage MakeTwoComponents() {
  return BinaryImage{.rows = 2, .cols = 4, .data = {1, 1, 0, 0, 1, 1, 0, 1}};
}

const std::array<TestCaseEntry, 7> kTestCaseTable = {{
    {.name = "empty", .make = MakeEmpty, .expected = {}},
    {.name = "single_pixel", .make = MakeSinglePixel, .expected = {Hull{Point{.x = 0, .y = 0}}}},
    {.name = "two_isolated",
     .make = MakeTwoIsolated,
     .expected = {Hull{Point{.x = 0, .y = 0}}, Hull{Point{.x = 1, .y = 1}}}},
    {.name = "vertical_line",
     .make = MakeVerticalLine,
     .expected = {Hull{Point{.x = 0, .y = 0}, Point{.x = 0, .y = 2}}}},
    {.name = "horizontal_line",
     .make = MakeHorizontalLine,
     .expected = {Hull{Point{.x = 0, .y = 0}, Point{.x = 2, .y = 0}}}},
    {.name = "filled_rect",
     .make = MakeFilledRect,
     .expected = {Hull{Point{.x = 0, .y = 0}, Point{.x = 2, .y = 0}, Point{.x = 2, .y = 1}, Point{.x = 0, .y = 1}}}},
    {.name = "two_components",
     .make = MakeTwoComponents,
     .expected = {Hull{Point{.x = 0, .y = 0}, Point{.x = 1, .y = 0}, Point{.x = 1, .y = 1}, Point{.x = 0, .y = 1}},
                  Hull{Point{.x = 3, .y = 1}}}},
}};

const TestCaseEntry *FindTestCase(const std::string &name) {
  for (const auto &e : kTestCaseTable) {
    if (name == e.name) {
      return &e;
    }
  }
  return nullptr;
}

}  // namespace

class KamalaginABinaryImageConvexHullFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return test_param;
  }

 protected:
  void SetUp() override {
    TestType name = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    const TestCaseEntry *entry = FindTestCase(name);
    ASSERT_NE(entry, nullptr);
    input_data_ = entry->make();
    expected_ = entry->expected;
    NormalizeHullList(expected_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    NormalizeHullList(output_data);
    return output_data == expected_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
  HullList expected_;
};

namespace {

TEST_P(KamalaginABinaryImageConvexHullFuncTests, RunCases) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 7> kTestNames = {"empty",           "single_pixel", "two_isolated",  "vertical_line",
                                            "horizontal_line", "filled_rect",  "two_components"};

const auto kTestTasksList = ppc::util::AddFuncTask<KamalaginABinaryImageConvexHullOMP, InType>(
    kTestNames, PPC_SETTINGS_kamalagin_a_binary_image_convex_hull_omp);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    KamalaginABinaryImageConvexHullFuncTests::PrintFuncTestName<KamalaginABinaryImageConvexHullFuncTests>;

INSTANTIATE_TEST_SUITE_P(BinaryImageConvexHullTests, KamalaginABinaryImageConvexHullFuncTests, kGtestValues,
                         kPerfTestName);

}  // namespace

}  // namespace kamalagin_a_binary_image_convex_hull_omp
