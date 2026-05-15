#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "peterson_r_graham_scan/common/include/common.hpp"
#include "peterson_r_graham_scan/seq/include/ops_seq.hpp"
#include "peterson_r_graham_scan/omp/include/ops_omp.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace peterson_r_graham_scan {

class PetersonGrahamScannerFuncTests : public ppc::util::BaseRunFuncTests<InputValue, OutputValue, TestParameters> {
 public:
  static std::string PrintTestParam(const TestParameters& test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestParameters params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
  }

  bool CheckTestOutputData(OutputValue& output_data) final {
    return input_data_ == output_data;
  }

  InputValue GetTestInputData() final {
    return input_data_;
  }

 private:
  InputValue input_data_ = 0;
};

namespace {

using Point = Point2D;

TEST_P(PetersonGrahamScannerFuncTests, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestParameters, 3> kTestCases = {std::make_tuple(3, "circle_3"), std::make_tuple(5, "circle_5"),
                                                  std::make_tuple(7, "circle_7")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<PetersonGrahamScannerSeq, InputValue>(kTestCases, PPC_SETTINGS_peterson_r_graham_scan),
    ppc::util::AddFuncTask<PetersonGrahamScannerOmp, InputValue>(kTestCases, PPC_SETTINGS_peterson_r_graham_scan));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kTestNameGenerator = PetersonGrahamScannerFuncTests::PrintFuncTestName<PetersonGrahamScannerFuncTests>;

INSTANTIATE_TEST_SUITE_P(DefaultTests, PetersonGrahamScannerFuncTests, kGtestValues, kTestNameGenerator);

}  // namespace

}  // namespace peterson_r_graham_scan
