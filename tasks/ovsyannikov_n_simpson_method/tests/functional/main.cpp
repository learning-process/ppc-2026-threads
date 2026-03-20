#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>

#include "ovsyannikov_n_simpson_method/common/include/common.hpp"
#include "ovsyannikov_n_simpson_method/omp/include/ops_omp.hpp"
#include "ovsyannikov_n_simpson_method/seq/include/ops_seq.hpp"
#include "ovsyannikov_n_simpson_method/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace ovsyannikov_n_simpson_method {

class OvsyannikovNRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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
    double dx = input_data_.bx - input_data_.ax;
    double dy = input_data_.by - input_data_.ay;
    double expected = dx * dy * (input_data_.ax + input_data_.bx + input_data_.ay + input_data_.by) / 2.0;

    return std::abs(output_data - expected) < 1e-4;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = {};
};

namespace {

TEST_P(OvsyannikovNRunFuncTestsThreads, SimpsonTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 6> kTestParam = {std::make_tuple(InType{0.0, 1.0, 0.0, 1.0, 10, 10}, "normal_10x10"),
                                            std::make_tuple(InType{0.0, 1.0, 0.0, 1.0, 50, 50}, "normal_50x50"),
                                            std::make_tuple(InType{1.0, 1.0, 0.0, 1.0, 10, 10}, "zero_width_x"),
                                            std::make_tuple(InType{0.0, 1.0, 5.0, 5.0, 10, 10}, "zero_width_y"),
                                            std::make_tuple(InType{-1.0, 0.0, 0.0, 1.0, 20, 20}, "negative_range_x"),
                                            std::make_tuple(InType{0.0, 2.0, 0.0, 2.0, 100, 50}, "large_range_2x2")};

const auto kTestTasksSEQ =
    ppc::util::AddFuncTask<OvsyannikovNSimpsonMethodSEQ, InType>(kTestParam, PPC_SETTINGS_ovsyannikov_n_simpson_method);
INSTANTIATE_TEST_SUITE_P(SimpsonTest_SEQ, OvsyannikovNRunFuncTestsThreads, ppc::util::ExpandToValues(kTestTasksSEQ),
                         OvsyannikovNRunFuncTestsThreads::PrintFuncTestName<OvsyannikovNRunFuncTestsThreads>);

const auto kTestTasksOMP =
    ppc::util::AddFuncTask<OvsyannikovNSimpsonMethodOMP, InType>(kTestParam, PPC_SETTINGS_ovsyannikov_n_simpson_method);
INSTANTIATE_TEST_SUITE_P(SimpsonTest_OMP, OvsyannikovNRunFuncTestsThreads, ppc::util::ExpandToValues(kTestTasksOMP),
                         OvsyannikovNRunFuncTestsThreads::PrintFuncTestName<OvsyannikovNRunFuncTestsThreads>);

const auto kTestTasksTBB =
    ppc::util::AddFuncTask<OvsyannikovNSimpsonMethodTBB, InType>(kTestParam, PPC_SETTINGS_ovsyannikov_n_simpson_method);
INSTANTIATE_TEST_SUITE_P(SimpsonTest_TBB, OvsyannikovNRunFuncTestsThreads, ppc::util::ExpandToValues(kTestTasksTBB),
                         OvsyannikovNRunFuncTestsThreads::PrintFuncTestName<OvsyannikovNRunFuncTestsThreads>);

}  // namespace
}  // namespace ovsyannikov_n_simpson_method
