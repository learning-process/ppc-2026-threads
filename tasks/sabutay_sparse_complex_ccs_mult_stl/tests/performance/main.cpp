#include <gtest/gtest.h>

#include <random>
#include <tuple>

#include "sabutay_sparse_complex_ccs_mult_stl/common/include/common.hpp"
#include "sabutay_sparse_complex_ccs_mult_stl/stl/include/ops_stl.hpp"
#include "util/include/perf_test_util.hpp"

namespace sabutay_sparse_complex_ccs_mult_stl {

namespace {

CCS CreateRandomSparseMatrix(int rows, int cols, double density = 0.1) {
  CCS matrix;
  matrix.m = rows;
  matrix.n = cols;
  matrix.col_ptr.assign(cols + 1, 0);

  std::random_device random_device;
  std::mt19937 random_generator(random_device());
  std::uniform_real_distribution<double> value_dist(-10.0, 10.0);
  std::uniform_int_distribution<int> row_dist(0, rows - 1);

  const int total_elements = static_cast<int>(rows * cols * density);
  for (int col = 0; col < cols; ++col) {
    const int elements_in_col = static_cast<int>((total_elements * (col + 1.0) / static_cast<double>(cols)) -
                                                 (total_elements * col / static_cast<double>(cols)));

    for (int i = 0; i < elements_in_col; ++i) {
      const int row = row_dist(random_generator);
      const double real_part = value_dist(random_generator);
      const double imag_part = value_dist(random_generator);

      matrix.row_ind.push_back(row);
      matrix.values.emplace_back(real_part, imag_part);
    }

    matrix.col_ptr[col + 1] = static_cast<int>(matrix.row_ind.size());
  }

  return matrix;
}

}  // namespace

class SabutayRunPerfTestsStl : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    matrix_a_ = CreateRandomSparseMatrix(100, 100, 0.05);
    matrix_b_ = CreateRandomSparseMatrix(100, 100, 0.05);
    input_data_ = std::make_tuple(matrix_a_, matrix_b_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const CCS &a = std::get<0>(input_data_);
    const CCS &b = std::get<1>(input_data_);
    return (output_data.m == a.m) && (output_data.n == b.n);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  CCS matrix_a_;
  CCS matrix_b_;
  InType input_data_;
};

TEST_P(SabutayRunPerfTestsStl, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, SabutaySparseComplexCcsMultSTL>(
    PPC_SETTINGS_sabutay_sparse_complex_ccs_mult_stl);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = SabutayRunPerfTestsStl::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, SabutayRunPerfTestsStl, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sabutay_sparse_complex_ccs_mult_stl
