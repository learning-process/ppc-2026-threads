#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <ranges>
#include <string>
#include <tuple>
#include <vector>

#include "util/include/func_test_util.hpp"
#include "yurkin_g_graham_scan/common/include/common.hpp"
#include "yurkin_g_graham_scan/seq/include/ops_seq.hpp"

namespace yurkin_g_graham_scan {

class YurkinGGrahamScanFuncTets : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    InType pts;
    pts.push_back({0.0, 0.0});
    pts.push_back({0.0, 1.0});
    pts.push_back({1.0, 0.0});
    pts.push_back({1.0, 1.0});
    pts.push_back({0.5, 0.5});  // interior
    pts.push_back({0.0, 0.0});  // duplicate
    pts.push_back({0.5, 0.0});
    pts.push_back({0.75, 0.0});
    input_data_ = pts;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    std::vector<Point> expected = {{0.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}, {1.0, 1.0}};
    if (output_data.size() != expected.size()) {
      return false;
    }

    auto contains = [](const std::vector<Point> &vec, const Point &p) {
      return std::ranges::any_of(vec, [&](const Point &q) { return q.x == p.x && q.y == p.y; });
    };

    return std::ranges::all_of(expected, [&](const Point &p) { return contains(output_data, p); });
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(YurkinGGrahamScanFuncTets, SquareWithInterior) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 1> kTestParam = {std::make_tuple(1, "square")};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<YurkinGGrahamScanSEQ, InType>(kTestParam, PPC_SETTINGS_yurkin_g_graham_scan));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = YurkinGGrahamScanFuncTets::PrintFuncTestName<YurkinGGrahamScanFuncTets>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, YurkinGGrahamScanFuncTets, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace yurkin_g_graham_scan
