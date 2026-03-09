#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>

#include "ovsyannikov_n_simpson_method_tbb/common/include/common.hpp"
#include "ovsyannikov_n_simpson_method_tbb/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace ovsyannikov_n_simpson_method_tbb {

class OvsyannikovNRunFuncTestsTBB : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    auto params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::abs(output_data - expected_val_) < 1e-4;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = {};
  double expected_val_ = 1.0;
};

namespace {

TEST_P(OvsyannikovNRunFuncTestsTBB, SimpsonTestTBB) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(InType{0.0, 1.0, 0.0, 1.0, 10, 10}, "steps_10"),
                                            std::make_tuple(InType{0.0, 1.0, 0.0, 1.0, 50, 50}, "steps_50"),
                                            std::make_tuple(InType{0.0, 1.0, 0.0, 1.0, 100, 100}, "steps_100")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<OvsyannikovNSimpsonMethodTBB, InType>(
    kTestParam, PPC_SETTINGS_ovsyannikov_n_simpson_method_omp));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = OvsyannikovNRunFuncTestsTBB::PrintFuncTestName<OvsyannikovNRunFuncTestsTBB>;
INSTANTIATE_TEST_SUITE_P(SimpsonTestTBB, OvsyannikovNRunFuncTestsTBB, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace ovsyannikov_n_simpson_method_tbb
