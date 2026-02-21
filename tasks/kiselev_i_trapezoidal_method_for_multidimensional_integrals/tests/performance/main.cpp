// #include <gtest/gtest.h>

// // #include "example_threads/all/include/ops_all.hpp"
// // #include "example_threads/common/include/common.hpp"
// // #include "example_threads/omp/include/ops_omp.hpp"
// // #include "example_threads/seq/include/ops_seq.hpp"
// // #include "example_threads/stl/include/ops_stl.hpp"
// // #include "example_threads/tbb/include/ops_tbb.hpp"

// #include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/common/include/common.hpp"
// #include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/seq/include/ops_seq.hpp"
// #include "util/include/perf_test_util.hpp"

// namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals {

// class KiselevIRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
//   const int kCount_ = 200;
//   InType input_data_{};

//   void SetUp() override {
//     input_data_ = kCount_;
//   }

//   bool CheckTestOutputData(OutType &output_data) final {
//     return input_data_ == output_data;
//   }

//   InType GetTestInputData() final {
//     return input_data_;
//   }
// };

// TEST_P(KiselevIRunPerfTestThreads, RunPerfModes) {
//   ExecuteTest(GetParam());
// }

// namespace {

// const auto kAllPerfTasks =
//     ppc::util::MakeAllPerfTasks<InType, KiselevITestTaskSEQ>(PPC_SETTINGS_example_threads);

// const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

// const auto kPerfTestName = KiselevIRunPerfTestThreads::CustomPerfTestName;

// INSTANTIATE_TEST_SUITE_P(RunModeTests, KiselevIRunPerfTestThreads, kGtestValues, kPerfTestName);

// }  // namespace

// }  // namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/common/include/common.hpp"
#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals {

class KiselevIRunPerfTestThreads
    : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  double expected_value_{};

  void SetUp() override {
    expected_value_ = 2.0 / 3.0;  // точный интеграл x^2 + y^2 на [0,1]^2
  }

  bool CheckTestOutputData(OutType& output_data) final {
    return std::abs(output_data - expected_value_) < 1e-2;
  }

  InType GetTestInputData() final {
    InType in;

    in.left_bounds = {0.0, 0.0};
    in.right_bounds = {1.0, 1.0};
    in.step_n_size = {200, 200};
    in.type_function = 0;   // x^2 + y^2
    in.epsilon = 0.0;       // фиксированная сетка

    return in;
  }
};

TEST_P(KiselevIRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KiselevITestTaskSEQ>(
        PPC_SETTINGS_example_threads);

const auto kGtestValues =
    ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName =
    KiselevIRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(
    RunModeTests,
    KiselevIRunPerfTestThreads,
    kGtestValues,
    kPerfTestName);

}  // namespace

}  // namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals