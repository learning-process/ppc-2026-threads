#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "perepelkin_i_convex_hull_graham_scan/common/include/common.hpp"
#include "perepelkin_i_convex_hull_graham_scan/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace perepelkin_i_convex_hull_graham_scan {

class PerepelkinIConvexHullGrahamScanFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<0>(test_param);
  }

 protected:
  void SetUp() override {
    const auto &[test_name, input_data, expected] =
        std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = input_data;
    expected_output_ = expected;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return expected_output_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_output_;
};

namespace {

TEST_P(PerepelkinIConvexHullGrahamScanFuncTests, ConvexHullGrahamScan) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 12> kTestParam = {
    std::make_tuple("empty_set", InType{}, OutType{}),
    std::make_tuple("single_point", InType{{{0.0, 0.0}}}, OutType{{{0.0, 0.0}}}),
    std::make_tuple("two_points", InType{{{1.0, 1.0}, {0.0, 0.0}}}, OutType{{{0.0, 0.0}, {1.0, 1.0}}}),
    std::make_tuple("collinear_points", InType{{{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {1.0, 1.0}}},
                    OutType{{{0.0, 0.0}, {2.0, 2.0}}}),
    std::make_tuple("rectangle", InType{{{1.0, 0.0}, {0.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}},
                    OutType{{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}}),
    std::make_tuple("duplicate_points", InType{{{0.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}}},
                    OutType{{{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}}}),
    std::make_tuple("same_angle_points", InType{{{0.0, 0.0}, {2.0, 2.0}, {1.0, 1.0}, {3.0, 3.0}, {0.0, 1.0}}},
                    OutType{{{0.0, 0.0}, {3.0, 3.0}, {0.0, 1.0}}}),
    std::make_tuple("complex_shape",
                    InType{{{0.0, 0.0}, {2.0, 1.0}, {1.0, -1.0}, {-1.0, -1.0}, {-2.0, 2.0}, {0.0, 2.0}, {1.0, 1.0}}},
                    OutType{{{-1.0, -1.0}, {1.0, -1.0}, {2.0, 1.0}, {0.0, 2.0}, {-2.0, 2.0}}}),
    std::make_tuple("decimal_coordinates", InType{{{0.1, 0.5}, {-0.2, 0.3}, {0.4, -0.6}, {1.1, 1.1}, {0.0, -1.0}}},
                    OutType{{{0.0, -1.0}, {0.4, -0.6}, {1.1, 1.1}, {0.1, 0.5}, {-0.2, 0.3}}}),
    std::make_tuple("irregular_positions",
                    InType{{{5.5, -3.2}, {-2.1, 4.8}, {3.3, 3.3}, {-4.4, -4.4}, {0.0, 2.2}, {-1.5, 0.0}}},
                    OutType{{{-4.4, -4.4}, {5.5, -3.2}, {3.3, 3.3}, {-2.1, 4.8}}}),
    std::make_tuple("huge_coordinates", InType{{{9e9, 9e9}, {-9e9, 9e9}, {9e9, -9e9}, {-9e9, -9e9}, {0.0, 0.0}}},
                    OutType{{{-9e9, -9e9}, {9e9, -9e9}, {9e9, 9e9}, {-9e9, 9e9}}}),
    std::make_tuple("many_inside",
                    InType{{{0.0, 0.0},
                            {1.0, 0.0},
                            {1.0, 1.0},
                            {0.0, 1.0},
                            {0.2, 0.2},
                            {0.3, 0.7},
                            {0.7, 0.3},
                            {0.6, 0.6},
                            {0.4, 0.5}}},
                    OutType{{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}})};

const auto kTestTasksList = ppc::util::AddFuncTask<PerepelkinIConvexHullGrahamScanSEQ, InType>(
    kTestParam, PPC_SETTINGS_perepelkin_i_convex_hull_graham_scan);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kFuncTestName =
    PerepelkinIConvexHullGrahamScanFuncTests::PrintFuncTestName<PerepelkinIConvexHullGrahamScanFuncTests>;

INSTANTIATE_TEST_SUITE_P(ConvexHullGrahamScanTests, PerepelkinIConvexHullGrahamScanFuncTests, kGtestValues,
                         kFuncTestName);

}  // namespace

}  // namespace perepelkin_i_convex_hull_graham_scan
