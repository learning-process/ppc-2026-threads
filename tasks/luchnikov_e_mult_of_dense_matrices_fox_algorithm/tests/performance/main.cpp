#include <gtest/gtest.h>

#include "luchnikov_e_mult_of_dense_matrices_fox_algorithm/common/include/common.hpp"
#include "luchnikov_e_mult_of_dense_matrices_fox_algorithm/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace luchnikov_e_mult_of_dense_matrices_fox_algorithm {

class LuchnikovEMultOfDenseMatrixFoxAlgoritmPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 200;
  InType input_data_{};

  void SetUp() override {
    input_data_ = kCount_;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data > 0;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(LuchnikovEMultOfDenseMatrixFoxAlgoritmPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {
const auto kAllPerfTasks =
    std::tuple_cat(ppc::util::MakeAllPerfTasks<InType, LuchnikovEMultOfDenseMatrixFoxAlgoritmSeq>(
                       PPC_SETTINGS_luchnikov_e_mult_of_dense_matrices_fox_algorithm),
                   ppc::util::MakeAllPerfTasks<InType, LuchnikovEMultOfDenseMatrixFoxAlgoritmOMP>(
                       PPC_SETTINGS_luchnikov_e_mult_of_dense_matrices_fox_algorithm),
                   ppc::util::MakeAllPerfTasks<InType, LuchnikovEMultOfDenseMatrixFoxAlgoritmSTL>(
                       PPC_SETTINGS_luchnikov_e_mult_of_dense_matrices_fox_algorithm),
                   ppc::util::MakeAllPerfTasks<InType, LuchnikovEMultOfDenseMatrixFoxAlgoritmTBB>(
                       PPC_SETTINGS_luchnikov_e_mult_of_dense_matrices_fox_algorithm));

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = LuchnikovEMultOfDenseMatrixFoxAlgoritmPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(AllTechPerfTests, LuchnikovEMultOfDenseMatrixFoxAlgoritmPerfTests, kGtestValues,
                         kPerfTestName);
}  // namespace
}  // namespace luchnikov_e_mult_of_dense_matrices_fox_algorithm
