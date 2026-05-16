#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>

#include "luchnikov_e_mult_of_dense_matrices_fox_algorithm_seq/common/include/common.hpp"
#include "luchnikov_e_mult_of_dense_matrices_fox_algorithm_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace luchnikov_e_mult_of_dense_matrices_fox_algorithm_seq {

class LuchnikovEMultOfDenseMatrixFoxAlgoritmPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  static bool CheckTestOutputData(OutType &output_data) final {
    return std::isfinite(static_cast<double>(output_data));
  }

  static InType GetTestInputData() final {
    static constexpr std::array<InType, 9> kBenchmarkSizes = {8, 12, 24, 48, 96, 192, 384, 512, 768};

    static thread_local size_t run_counter = 0;
    InType size = kBenchmarkSizes.at(run_counter % kBenchmarkSizes.size());
    ++run_counter;
    return size;
  }
};

static TestP(LuchnikovEMultOfDenseMatrixFoxAlgoritmPerfTestThreads /*unused*/, BenchmarkFoxAlgorithm /*unused*/) {
  ExecuteTest(GetParam());
}

namespace {

const auto kPerfTaskList =
    std::tuple_cat(ppc::util::MakeAllPerfTasks<InType, LuchnikovEMultOfDenseMatrixFoxAlgoritmSeq>(
        PPC_SETTINGS_luchnikov_e_mult_of_dense_matrices_fox_algorithm));

const auto kGtestParams = ppc::util::TupleToGTestValues(kPerfTaskList);
const auto kTestPrinter = LuchnikovEMultOfDenseMatrixFoxAlgoritmPerfTestThreads::CustomPerfTestName;

InstantiateTestSuiteP(FoxAlgorithmBenchmark, LuchnikovEMultOfDenseMatrixFoxAlgoritmPerfTestThreads, kGtestParams,
                      kTestPrinter);

}  // namespace
}  // namespace luchnikov_e_mult_of_dense_matrices_fox_algorithm_seq
