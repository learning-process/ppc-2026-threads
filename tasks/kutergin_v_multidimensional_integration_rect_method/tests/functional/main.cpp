#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "../../all/include/rect_method_all.hpp"
#include "../../common/include/common.hpp"
#include "../../omp/include/rect_method_omp.hpp"
#include "../../seq/include/rect_method_seq.hpp"
#include "../../stl/include/rect_method_stl.hpp"
#include "../../tbb/include/rect_method_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kutergin_v_multidimensional_integration_rect_method {

class RectMethodFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    const auto &params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_output_ = std::get<1>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    int rank = 0;
    if (ppc::util::IsUnderMpirun()) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }

    if (rank != 0) {
      return true;
    }

    const auto &params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    std::string test_name = std::get<2>(params);

    std::cout << "[ DEBUG ] Test: " << test_name << '\n';
    std::cout << "[ DEBUG ] Expected: " << expected_output_ << '\n';
    std::cout << "[ DEBUG ] Actual:   " << output_data << '\n';
    std::cout << "[ DEBUG ] Diff:     " << std::abs(output_data - expected_output_) << '\n';
    return std::abs(output_data - expected_output_) < 1e-4;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
  OutType expected_output_{0.0};
};

namespace {

TEST_P(RectMethodFuncTests, IntegrationTest) {
  ExecuteTest(GetParam());
}

IntegrationData test_1d = {
    .limits = {{0.0, 1.0}}, .n_steps = {100}, .func = [](const std::vector<double> &x) { return x[0]; }};

IntegrationData test_2d = {.limits = {{0.0, 1.0}, {0.0, 1.0}},
                           .n_steps = {50, 50},
                           .func = [](const std::vector<double> &x) { return x[0] + x[1]; }};

IntegrationData test_3d = {.limits = {{0.0, 2.0}, {0.0, 3.0}, {0.0, 1.0}},
                           .n_steps = {20, 20, 20},
                           .func = [](const std::vector<double> & /*unused*/) { return 1.0; }};

IntegrationData test_4d = {.limits = {{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}},
                           .n_steps = {5, 5, 5, 5},
                           .func = [](const std::vector<double> &x) { return x[0] + x[1] + x[2] + x[3]; }};

IntegrationData test_5d = {.limits = {{0.0, 2.0}, {0.0, 3.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}},
                           .n_steps = {3, 3, 3, 3, 3},
                           .func = [](const std::vector<double> &x) { return x[0] + x[1] + x[2] + x[3] + x[4]; }};

const std::array<TestType, 5> kTestCases = {
    std::make_tuple(test_1d, 0.5, "1D_Linear"), std::make_tuple(test_2d, 1.0, "2D_Linear"),
    std::make_tuple(test_3d, 6.0, "3D_Constant"), std::make_tuple(test_4d, 2.0, "4D_Sum_Small_Grid"),
    std::make_tuple(test_5d, 24.0, "5D_Sum_Small_Grid")};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<RectMethodSequential, InType>(
                       kTestCases, PPC_SETTINGS_kutergin_v_multidimensional_integration_rect_method),
                   ppc::util::AddFuncTask<RectMethodOMP, InType>(
                       kTestCases, PPC_SETTINGS_kutergin_v_multidimensional_integration_rect_method),
                   ppc::util::AddFuncTask<RectMethodTBB, InType>(
                       kTestCases, PPC_SETTINGS_kutergin_v_multidimensional_integration_rect_method),
                   ppc::util::AddFuncTask<RectMethodSTL, InType>(
                       kTestCases, PPC_SETTINGS_kutergin_v_multidimensional_integration_rect_method),
                   ppc::util::AddFuncTask<RectMethodALL, InType>(
                       kTestCases, PPC_SETTINGS_kutergin_v_multidimensional_integration_rect_method));

const auto kGTestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = RectMethodFuncTests::PrintFuncTestName<RectMethodFuncTests>;

INSTANTIATE_TEST_SUITE_P(MultidimensionalIntegration, RectMethodFuncTests, kGTestValues, kTestName);

}  // namespace

}  // namespace kutergin_v_multidimensional_integration_rect_method
