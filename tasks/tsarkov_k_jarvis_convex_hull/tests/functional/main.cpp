#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "tsarkov_k_jarvis_convex_hull/common/include/common.hpp"
#include "tsarkov_k_jarvis_convex_hull/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace tsarkov_k_jarvis_convex_hull {

class TsarkovKRunFuncTestsSEQ : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    const int test_case_id =
        std::get<0>(std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam()));

    if (test_case_id == 0) {
      input_data_ = {Point{.x = 0, .y = 0}};
      expected_output_ = {Point{.x = 0, .y = 0}};
    } else if (test_case_id == 1) {
      input_data_ = {Point{.x = 0, .y = 0}, Point{.x = 2, .y = 0}};
      expected_output_ = {Point{.x = 0, .y = 0}, Point{.x = 2, .y = 0}};
    } else if (test_case_id == 2) {
      input_data_ = {Point{.x = 0, .y = 0}, Point{.x = 2, .y = 0}, Point{.x = 1, .y = 2}};
      expected_output_ = {Point{.x = 0, .y = 0}, Point{.x = 2, .y = 0}, Point{.x = 1, .y = 2}};
    } else if (test_case_id == 3) {
      input_data_ = {Point{.x = 0, .y = 0}, Point{.x = 2, .y = 0}, Point{.x = 2, .y = 2}, Point{.x = 0, .y = 2},
                     Point{.x = 1, .y = 1}};
      expected_output_ = {Point{.x = 0, .y = 0}, Point{.x = 2, .y = 0}, Point{.x = 2, .y = 2}, Point{.x = 0, .y = 2}};
    } else if (test_case_id == 4) {
      input_data_ = {Point{.x = 0, .y = 0}, Point{.x = 1, .y = 1}, Point{.x = 2, .y = 2}, Point{.x = 3, .y = 3},
                     Point{.x = 4, .y = 4}};
      expected_output_ = {Point{.x = 0, .y = 0}, Point{.x = 4, .y = 4}};
    } else if (test_case_id == 5) {
      input_data_ = {Point{.x = 0, .y = 0}, Point{.x = 2, .y = 0}, Point{.x = 2, .y = 2}, Point{.x = 0, .y = 2},
                     Point{.x = 0, .y = 0}, Point{.x = 2, .y = 0}, Point{.x = 1, .y = 1}, Point{.x = 1, .y = 1}};
      expected_output_ = {Point{.x = 0, .y = 0}, Point{.x = 2, .y = 0}, Point{.x = 2, .y = 2}, Point{.x = 0, .y = 2}};
    } else if (test_case_id == 6) {
      input_data_ = {Point{.x = -2, .y = 0}, Point{.x = -1, .y = -1}, Point{.x = 0, .y = -2}, Point{.x = 2, .y = 0},
                     Point{.x = 0, .y = 2},  Point{.x = -1, .y = 1},  Point{.x = 0, .y = 0}};
      expected_output_ = {Point{.x = -2, .y = 0}, Point{.x = 0, .y = -2}, Point{.x = 2, .y = 0}, Point{.x = 0, .y = 2}};
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == expected_output_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_output_;
};

namespace {

TEST_P(TsarkovKRunFuncTestsSEQ, RunJarvisConvexHull) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 7> kTestParam = {
    std::make_tuple(0, "single_point"),       std::make_tuple(1, "two_points"),
    std::make_tuple(2, "triangle"),           std::make_tuple(3, "square_with_inner_point"),
    std::make_tuple(4, "collinear_points"),   std::make_tuple(5, "duplicate_points"),
    std::make_tuple(6, "diamond_with_inside")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<TsarkovKJarvisConvexHullSEQ, InType>(kTestParam, PPC_SETTINGS_tsarkov_k_jarvis_convex_hull));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = TsarkovKRunFuncTestsSEQ::PrintFuncTestName<TsarkovKRunFuncTestsSEQ>;

INSTANTIATE_TEST_SUITE_P(JarvisHullTests, TsarkovKRunFuncTestsSEQ, kGtestValues, kTestName);

}  // namespace

}  // namespace tsarkov_k_jarvis_convex_hull
