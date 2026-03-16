#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "dergachev_a_graham_scan/common/include/common.hpp"
#include "dergachev_a_graham_scan/omp/include/ops_omp.hpp"
#include "dergachev_a_graham_scan/seq/include/ops_seq.hpp"
#include "dergachev_a_graham_scan/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace dergachev_a_graham_scan {

class DergachevAGrahamScanFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
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

TEST_P(DergachevAGrahamScanFuncTestsThreads, GrahamScan) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 13> kTestParam = {
    std::make_tuple(0, "empty"),         std::make_tuple(1, "single_point"), std::make_tuple(2, "two_points"),
    std::make_tuple(3, "circle_3"),      std::make_tuple(4, "circle_4"),     std::make_tuple(5, "circle_5"),
    std::make_tuple(6, "circle_6"),      std::make_tuple(7, "circle_7"),     std::make_tuple(10, "circle_10"),
    std::make_tuple(15, "circle_15"),    std::make_tuple(50, "circle_50"),   std::make_tuple(100, "circle_100"),
    std::make_tuple(1000, "circle_1000")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<DergachevAGrahamScanOMP, InType>(kTestParam, PPC_SETTINGS_dergachev_a_graham_scan),
    ppc::util::AddFuncTask<DergachevAGrahamScanSEQ, InType>(kTestParam, PPC_SETTINGS_dergachev_a_graham_scan),
    ppc::util::AddFuncTask<DergachevAGrahamScanTBB, InType>(kTestParam, PPC_SETTINGS_dergachev_a_graham_scan));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    DergachevAGrahamScanFuncTestsThreads::PrintFuncTestName<DergachevAGrahamScanFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(GrahamScanTests, DergachevAGrahamScanFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace dergachev_a_graham_scan
