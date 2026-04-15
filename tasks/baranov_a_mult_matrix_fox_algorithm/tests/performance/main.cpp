#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "baranov_a_mult_matrix_fox_algorithm/all/include/ops_all.hpp"
#include "baranov_a_mult_matrix_fox_algorithm/common/include/common.hpp"
#include "baranov_a_mult_matrix_fox_algorithm/omp/include/ops_omp.hpp"
#include "baranov_a_mult_matrix_fox_algorithm/seq/include/ops_seq.hpp"
#include "baranov_a_mult_matrix_fox_algorithm/stl/include/ops_stl.hpp"
#include "baranov_a_mult_matrix_fox_algorithm/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace baranov_a_mult_matrix_fox_algorithm_test {

class BaranovAPerfTests : public ppc::util::BaseRunPerfTests<baranov_a_mult_matrix_fox_algorithm::InType,
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

TEST_P(BaranovAPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = std::tuple_cat(
    ppc::util::MakeAllPerfTasks<baranov_a_mult_matrix_fox_algorithm::InType,
                                baranov_a_mult_matrix_fox_algorithm_seq::BaranovAMultMatrixFoxAlgorithmSEQ>(
        PPC_SETTINGS_baranov_a_mult_matrix_fox_algorithm),
    ppc::util::MakeAllPerfTasks<baranov_a_mult_matrix_fox_algorithm::InType,
                                baranov_a_mult_matrix_fox_algorithm_omp::BaranovAMultMatrixFoxAlgorithmOMP>(
        PPC_SETTINGS_baranov_a_mult_matrix_fox_algorithm),
    ppc::util::MakeAllPerfTasks<baranov_a_mult_matrix_fox_algorithm::InType,
                                baranov_a_mult_matrix_fox_algorithm_tbb::BaranovAMultMatrixFoxAlgorithmTBB>(
        PPC_SETTINGS_baranov_a_mult_matrix_fox_algorithm),
    ppc::util::MakeAllPerfTasks<baranov_a_mult_matrix_fox_algorithm::InType,
                                baranov_a_mult_matrix_fox_algorithm_stl::BaranovAMultMatrixFoxAlgorithmSTL>(
        PPC_SETTINGS_baranov_a_mult_matrix_fox_algorithm),
    ppc::util::MakeAllPerfTasks<baranov_a_mult_matrix_fox_algorithm::InType,
                                baranov_a_mult_matrix_fox_algorithm_all::BaranovAMultMatrixFoxAlgorithmALL>(
        PPC_SETTINGS_baranov_a_mult_matrix_fox_algorithm));

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = BaranovAPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, BaranovAPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace baranov_a_mult_matrix_fox_algorithm_test
