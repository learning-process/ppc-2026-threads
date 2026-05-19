#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "egorova_l_binary_convex_hull/common/include/common.hpp"
#include "egorova_l_binary_convex_hull/omp/include/ops_omp.hpp"
#include "egorova_l_binary_convex_hull/seq/include/ops_seq.hpp"
#include "egorova_l_binary_convex_hull/stl/include/ops_stl.hpp"
#include "egorova_l_binary_convex_hull/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace egorova_l_binary_convex_hull {

using TestType = std::tuple<InType, std::vector<std::vector<Point>>, std::string>;

class EgorovaLFuncTest : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_result_ = std::get<1>(params);
  }

  static bool HullsEqual(const std::vector<Point> &lhs, const std::vector<Point> &rhs) {
    if (lhs.size() != rhs.size()) {
      return false;
    }
    for (size_t j = 0; j < lhs.size(); ++j) {
      if (lhs[j].x != rhs[j].x || lhs[j].y != rhs[j].y) {
        return false;
      }
    }
    return true;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_result_.size()) {
      return false;
    }
    auto sorted_output = output_data;
    auto sorted_expected = expected_result_;

    auto sort_hulls = [](std::vector<std::vector<Point>> &hulls) {
      for (auto &h : hulls) {
        std::ranges::sort(h, [](const Point &a, const Point &b) { return std::tie(a.x, a.y) < std::tie(b.x, b.y); });
      }
      std::ranges::sort(hulls, [](const std::vector<Point> &a, const std::vector<Point> &b) {
        if (a.empty() || b.empty()) {
          return a.size() < b.size();
        }
        return std::tie(a[0].x, a[0].y) < std::tie(b[0].x, b[0].y);
      });
    };

    sort_hulls(sorted_output);
    sort_hulls(sorted_expected);

    for (size_t i = 0; i < sorted_output.size(); ++i) {
      if (!HullsEqual(sorted_output[i], sorted_expected[i])) {
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
  img.data.assign(static_cast<size_t>(width) * height, 0);
  return img;
}

void DrawRect(InType &img, int x1, int y1, int x2, int y2) {
  for (int py = y1; py <= y2; ++py) {
    for (int px = x1; px <= x2; ++px) {
      if (px >= 0 && px < img.width && py >= 0 && py < img.height) {
        img.data[(static_cast<size_t>(py) * img.width) + px] = 255;
      }
    }
  }
}

const std::array<TestType, 2> kParams = {
    {std::make_tuple(
         []() {
  auto img = CreateEmptyImage(10, 10);
  DrawRect(img, 1, 1, 3, 3);
  return img;
}(), std::vector<std::vector<Point>>{{{1, 1}, {3, 1}, {3, 3}, {1, 3}}}, "single_sq"),
     std::make_tuple(CreateEmptyImage(10, 10), std::vector<std::vector<Point>>{}, "empty")}};

const auto kTaskList = std::tuple_cat(
    ppc::util::AddFuncTask<BinaryConvexHullSEQ, InType>(kParams, PPC_SETTINGS_egorova_l_binary_convex_hull),
    ppc::util::AddFuncTask<BinaryConvexHullOMP, InType>(kParams, PPC_SETTINGS_egorova_l_binary_convex_hull),
    ppc::util::AddFuncTask<BinaryConvexHullTBB, InType>(kParams, PPC_SETTINGS_egorova_l_binary_convex_hull),
    ppc::util::AddFuncTask<BinaryConvexHullSTL, InType>(kParams, PPC_SETTINGS_egorova_l_binary_convex_hull));

TEST_P(EgorovaLFuncTest, RunFunctionalTests) {
  ExecuteTest(GetParam());
}

INSTANTIATE_TEST_SUITE_P(BinaryConvexHullTests, EgorovaLFuncTest, ppc::util::ExpandToValues(kTaskList),
                         EgorovaLFuncTest::PrintFuncTestName<EgorovaLFuncTest>);

}  // namespace
}  // namespace egorova_l_binary_convex_hull
