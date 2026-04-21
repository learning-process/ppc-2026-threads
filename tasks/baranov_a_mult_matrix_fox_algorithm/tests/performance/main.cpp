#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "baranov_a_mult_matrix_fox_algorithm/common/include/common.hpp"
#include "baranov_a_mult_matrix_fox_algorithm/omp/include/ops_omp.hpp"
#include "baranov_a_mult_matrix_fox_algorithm/seq/include/ops_seq.hpp"
#include "baranov_a_mult_matrix_fox_algorithm/stl/include/ops_stl.hpp"
#include "baranov_a_mult_matrix_fox_algorithm/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace baranov_a_mult_matrix_fox_algorithm_test {

template <typename TaskType>
class BaranovAMultMatrixFoxAlgorithmPerfTests
    : public ppc::util::BaseRunPerfTests<baranov_a_mult_matrix_fox_algorithm::InType,
                                         baranov_a_mult_matrix_fox_algorithm::OutType> {
  void SetUp() override {
    size_t n = 512;

    size_t size = n * n;
    std::vector<double> a(size, 1.5);
    std::vector<double> b(size, 2.0);

    input_data_ = std::make_tuple(n, a, b);

    double expected_value = 3.0 * static_cast<double>(n);
    expected_output_.resize(size, expected_value);
  }

  bool CheckTestOutputData(baranov_a_mult_matrix_fox_algorithm::OutType &output_data) final {
    if (expected_output_.size() != output_data.size()) {
      return false;
    }

    double epsilon = 1e-8;
    for (size_t i = 0; i < expected_output_.size(); ++i) {
      if (std::abs(expected_output_[i] - output_data[i]) > epsilon) {
        return false;
      }
    }
    return true;
  }

  baranov_a_mult_matrix_fox_algorithm::InType GetTestInputData() final {
    return input_data_;
  }

 private:
  baranov_a_mult_matrix_fox_algorithm::InType input_data_;
  baranov_a_mult_matrix_fox_algorithm::OutType expected_output_;
};

using BaranovASEQPerfTests =
    BaranovAMultMatrixFoxAlgorithmPerfTests<baranov_a_mult_matrix_fox_algorithm_seq::BaranovAMultMatrixFoxAlgorithmSEQ>;
using BaranovAOMPPerfTests =
    BaranovAMultMatrixFoxAlgorithmPerfTests<baranov_a_mult_matrix_fox_algorithm_omp::BaranovAMultMatrixFoxAlgorithmOMP>;
using BaranovATBBPerfTests =
    BaranovAMultMatrixFoxAlgorithmPerfTests<baranov_a_mult_matrix_fox_algorithm_tbb::BaranovAMultMatrixFoxAlgorithmTBB>;
using BaranovASTLPerfTests =
    BaranovAMultMatrixFoxAlgorithmPerfTests<baranov_a_mult_matrix_fox_algorithm_stl::BaranovAMultMatrixFoxAlgorithmSTL>;

TEST_P(BaranovASEQPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

TEST_P(BaranovAOMPPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

TEST_P(BaranovATBBPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

TEST_P(BaranovASTLPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasksSeq =
    ppc::util::MakeAllPerfTasks<baranov_a_mult_matrix_fox_algorithm::InType,
                                baranov_a_mult_matrix_fox_algorithm_seq::BaranovAMultMatrixFoxAlgorithmSEQ>(
        PPC_SETTINGS_baranov_a_mult_matrix_fox_algorithm);

const auto kAllPerfTasksOmp =
    ppc::util::MakeAllPerfTasks<baranov_a_mult_matrix_fox_algorithm::InType,
                                baranov_a_mult_matrix_fox_algorithm_omp::BaranovAMultMatrixFoxAlgorithmOMP>(
        PPC_SETTINGS_baranov_a_mult_matrix_fox_algorithm);

const auto kAllPerfTasksTbb =
    ppc::util::MakeAllPerfTasks<baranov_a_mult_matrix_fox_algorithm::InType,
                                baranov_a_mult_matrix_fox_algorithm_tbb::BaranovAMultMatrixFoxAlgorithmTBB>(
        PPC_SETTINGS_baranov_a_mult_matrix_fox_algorithm);

const auto kAllPerfTasksStl =
    ppc::util::MakeAllPerfTasks<baranov_a_mult_matrix_fox_algorithm::InType,
                                baranov_a_mult_matrix_fox_algorithm_stl::BaranovAMultMatrixFoxAlgorithmSTL>(
        PPC_SETTINGS_baranov_a_mult_matrix_fox_algorithm);

const auto kGtestValuesSeq = ppc::util::TupleToGTestValues(kAllPerfTasksSeq);
const auto kGtestValuesOmp = ppc::util::TupleToGTestValues(kAllPerfTasksOmp);
const auto kGtestValuesTbb = ppc::util::TupleToGTestValues(kAllPerfTasksTbb);
const auto kGtestValuesStl = ppc::util::TupleToGTestValues(kAllPerfTasksStl);

const auto kPerfTestNameSeq = BaranovASEQPerfTests::CustomPerfTestName;
const auto kPerfTestNameOmp = BaranovAOMPPerfTests::CustomPerfTestName;
const auto kPerfTestNameTbb = BaranovATBBPerfTests::CustomPerfTestName;
const auto kPerfTestNameStl = BaranovASTLPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerfTestsSeq, BaranovASEQPerfTests, kGtestValuesSeq, kPerfTestNameSeq);
INSTANTIATE_TEST_SUITE_P(PerfTestsOmp, BaranovAOMPPerfTests, kGtestValuesOmp, kPerfTestNameOmp);
INSTANTIATE_TEST_SUITE_P(PerfTestsTbb, BaranovATBBPerfTests, kGtestValuesTbb, kPerfTestNameTbb);
INSTANTIATE_TEST_SUITE_P(PerfTestsStl, BaranovASTLPerfTests, kGtestValuesStl, kPerfTestNameStl);

}  // namespace

}  // namespace baranov_a_mult_matrix_fox_algorithm_test
