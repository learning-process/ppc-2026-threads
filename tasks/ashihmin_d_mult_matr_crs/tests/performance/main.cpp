#pragma once

#include <gtest/gtest.h>

#include <algorithm>
#include <tuple>
#include <vector>

#include "ashihmin_d_mult_matr_crs/common/include/common.hpp"
#include "ashihmin_d_mult_matr_crs/seq/include/ops_seq.hpp"
#include "performance/include/performance.hpp"
#include "util/include/perf_test_util.hpp"

namespace ashihmin_d_mult_matr_crs {

namespace {

CRSMatrix GenerateBandMatrix(std::size_t matrix_size, std::size_t bandwidth, double fill_value) {
  CRSMatrix result_matrix;
  result_matrix.rows = matrix_size;
  result_matrix.cols = matrix_size;
  result_matrix.row_ptr.resize(matrix_size + 1, 0);

  for (std::size_t row_index = 0; row_index < matrix_size; ++row_index) {
    const std::size_t begin_index = (row_index > bandwidth ? row_index - bandwidth : 0);
    const std::size_t end_index = std::min(matrix_size - 1, row_index + bandwidth);

    for (std::size_t column_index = begin_index; column_index <= end_index; ++column_index) {
      result_matrix.values.push_back(fill_value);
      result_matrix.col_index.push_back(static_cast<int>(column_index));
    }

    result_matrix.row_ptr[row_index + 1] = result_matrix.values.size();
  }

  return result_matrix;
}

}  // namespace

class AshihminDMultMatrCrsPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  static constexpr std::size_t kMatrixSize = 40000;
  static constexpr std::size_t kBandwidth = 30;

  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) override {
    ppc::util::BaseRunPerfTests<InType, OutType>::SetPerfAttributes(perf_attrs);
    perf_attrs.num_running = 1;
  }

  void SetUp() override {
    CRSMatrix matrix_a = GenerateBandMatrix(kMatrixSize, kBandwidth, 2.0);
    CRSMatrix matrix_b = GenerateBandMatrix(kMatrixSize, kBandwidth, 3.0);
    input_data_ = std::make_tuple(matrix_a, matrix_b);

    expected_.rows = kMatrixSize;
    expected_.cols = kMatrixSize;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.rows == expected_.rows && output_data.cols == expected_.cols && !output_data.values.empty();
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  CRSMatrix expected_;
};

TEST_P(AshihminDMultMatrCrsPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, AshihminDMultMatrCrsSEQ>(PPC_SETTINGS_ashihmin_d_mult_matr_crs);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = AshihminDMultMatrCrsPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(AshihminSparseCRSPerfTests, AshihminDMultMatrCrsPerfTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace ashihmin_d_mult_matr_crs
