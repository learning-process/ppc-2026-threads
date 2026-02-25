#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "perepelkin_i_convex_hull_graham_scan/common/include/common.hpp"
// #include "perepelkin_i_convex_hull_graham_scan/mpi/include/ops_mpi.hpp"
#include "perepelkin_i_convex_hull_graham_scan/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace perepelkin_i_convex_hull_graham_scan {

class PerepelkinIConvexHullGrahamScanFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);   // update
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    // input_data_ = std::get<0>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (input_data_ == output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = 0;
};

namespace {

TEST_P(PerepelkinIConvexHullGrahamScanFuncTests, ConvexHullGrahamScan) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {
  std::make_tuple(3, "3"), 
  std::make_tuple(5, "5"), 
  std::make_tuple(7, "7")
};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<PerepelkinIConvexHullGrahamScanSEQ, InType>(kTestParam, PPC_SETTINGS_perepelkin_i_convex_hull_graham_scan));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kFuncTestName = PerepelkinIConvexHullGrahamScanFuncTests::PrintFuncTestName<PerepelkinIConvexHullGrahamScanFuncTests>;

INSTANTIATE_TEST_SUITE_P(ConvexHullGrahamScanTests, PerepelkinIConvexHullGrahamScanFuncTests, kGtestValues, kFuncTestName);

}  // namespace

}  // namespace perepelkin_i_convex_hull_graham_scan
