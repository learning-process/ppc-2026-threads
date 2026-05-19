#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <cstddef>
#include <tuple>
#include "util/include/perf_test_util.hpp"
#include "luchnikov_e_mult_of_dense_matrices_fox_algorithm_tbb/common/include/common.hpp"
#include "luchnikov_e_mult_of_dense_matrices_fox_algorithm_tbb/tbb/include/ops_tbb.hpp"

namespace luchnikov_e_mult_of_dense_matrices_fox_algorithm_tbb {

class LuchnikovEMultOfDenseMatrixFoxAlgoritmPerfTestTBB
    : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  bool CheckTestOutputData(OutType &output_data) final {
    return std::isfinite(static_cast<double>(output_data));
  }

  InType GetTestInputData() final {
    static constexpr std::array<InType, 9> kBenchmarkSizes = {
        8, 12, 24, 48, 96, 192, 384, 512, 768
    };
    
    static thread_local size_t run_counter = 0;
    InType size = kBenchmarkSizes.at(run_counter % kBenchmarkSizes.size());
    ++run_counter;
    return size;
  }
};

TEST_P(LuchnikovEMultOfDenseMatrixFoxAlgoritmPerfTestTBB, BenchmarkFoxAlgorithm) {
  ExecuteTest(GetParam());
}

namespace {

const auto kPerfTaskList = std::tuple_cat(
    ppc::util::MakeAllPerfTasks<InType, LuchnikovEMultOfDenseMatrixFoxAlgoritmTBB>(
        PPC_SETTINGS_luchnikov_e_mult_of_dense_matrices_fox_algorithm_tbb));

const auto kGtestParams = ppc::util::TupleToGTestValues(kPerfTaskList);
const auto kTestPrinter = LuchnikovEMultOfDenseMatrixFoxAlgoritmPerfTestTBB::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(
    TBBFoxAlgorithmBenchmark,
    LuchnikovEMultOfDenseMatrixFoxAlgoritmPerfTestTBB,
    kGtestParams,
    kTestPrinter);

}  // namespace
}  // namespace luchnikov_e_mult_of_dense_matrices_fox_algorithm_tbb