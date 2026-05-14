#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "dorogin_v_bin_img_conv_hull_stl/common/include/common.hpp"
#include "dorogin_v_bin_img_conv_hull_stl/stl/include/ops_stl.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace nesterov_a_test_task_threads {

class NesterovARunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }

  static InType MakeImage(const int w, const int h, const std::vector<std::uint8_t> &data) {
    return InType{.width = w, .height = h, .data = data};
  }

 protected:
  static std::vector<Point> NormalizeHull(std::vector<Point> hull) {
    std::ranges::sort(hull,
                      [](const Point &a, const Point &b) { return (a.x < b.x) || ((a.x == b.x) && (a.y < b.y)); });
    return hull;
  }

  static OutType NormalizeAll(OutType hulls) {
    for (auto &h : hulls) {
      h = NormalizeHull(std::move(h));
    }
    std::ranges::sort(hulls, [](const std::vector<Point> &a, const std::vector<Point> &b) {
      if (a.empty() || b.empty()) {
        return a.size() < b.size();
      }
      if (a.front().x != b.front().x) {
        return a.front().x < b.front().x;
      }
      return a.front().y < b.front().y;
    });
    return hulls;
  }

  void SetUp() override {
    const auto params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_tag_ = std::get<1>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    OutType expected{};
    if (expected_tag_ == "single_pixel") {
      expected = {{{Point{.x = 1, .y = 1}}}};
    } else if (expected_tag_ == "filled_square") {
      expected = {{{Point{.x = 1, .y = 1}, Point{.x = 3, .y = 1}, Point{.x = 3, .y = 3}, Point{.x = 1, .y = 3}}}};
    } else if (expected_tag_ == "two_components") {
      expected = {{{Point{.x = 0, .y = 0}, Point{.x = 1, .y = 0}, Point{.x = 1, .y = 1}, Point{.x = 0, .y = 1}}},
                  {{Point{.x = 4, .y = 3}, Point{.x = 5, .y = 3}, Point{.x = 5, .y = 4}, Point{.x = 4, .y = 4}}}};
    } else {
      return false;
    }
    return NormalizeAll(std::move(output_data)) == NormalizeAll(std::move(expected));
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
  std::string expected_tag_;
};

namespace {

TEST_P(NesterovARunFuncTestsThreads, ConvexHullForBinaryComponentsSTL) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {
    std::make_tuple(NesterovARunFuncTestsThreads::MakeImage(3, 3, {0, 0, 0, 0, 1, 0, 0, 0, 0}), "single_pixel"),
    std::make_tuple(NesterovARunFuncTestsThreads::MakeImage(
                        5, 5, {0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0}),
                    "filled_square"),
    std::make_tuple(NesterovARunFuncTestsThreads::MakeImage(6, 5, {1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                                                                   0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1}),
                    "two_components")};

const auto kStlTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<NesterovATestTaskSTL, InType>(kTestParam, PPC_SETTINGS_dorogin_v_bin_img_conv_hull_stl));
const auto kStlValues = ppc::util::ExpandToValues(kStlTasksList);

const auto kNameFn = NesterovARunFuncTestsThreads::PrintFuncTestName<NesterovARunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(BinaryImageConvexHullStlTests, NesterovARunFuncTestsThreads, kStlValues, kNameFn);

}  // namespace

}  // namespace nesterov_a_test_task_threads
