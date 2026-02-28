#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <vector>

#include "ashihmin_d_mult_matr_crs/common/include/common.hpp"
#include "ashihmin_d_mult_matr_crs/seq/include/ops_seq.hpp"
#include "performance/include/performance.hpp"
#include "util/include/perf_test_util.hpp"

namespace ashihmin_d_mult_matr_crs {

namespace {

CRSMatrix GenerateBandMatrix(std::size_t matrix_size, std::size_t bandwidth, double fill_value) {
  CRSMatrix matrix_result;
  matrix_result.rows = matrix_size;
  matrix_result.cols = matrix_size;
  matrix_result.row_ptr.resize(matrix_size + 1, 0);

  for (std::size_t row_index = 0; row_index < matrix_size; ++row_index) {
    std::size_t begin_col = (row_index > bandwidth ? row_index - bandwidth : 0);
    std::size_t end_col = std::min(matrix_size - 1, row_index + bandwidth);

    for (std::size_t col_index = begin_col; col_index <= end_col; ++col_index) {
      matrix_result.values.push_back(fill_value);
      matrix_result.col_index.push_back(static_cast<int>(col_index));
    }

    matrix_result.row_ptr[row_index + 1] = static_cast<int>(matrix_result.values.size());
  }

  return matrix_result;
}

}  // namespace

class AshihminDMultMatrCrsPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  static constexpr std::size_t kMatrixSize = 40000;
  static constexpr std::size_t kMatrixBandwidth = 30;

  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attributes) override {
    ppc::util::BaseRunPerfTests<InType, OutType>::SetPerfAttributes(perf_attributes);
    perf_attributes.num_running = 1;  // одно тяжелое прогон, чтобы метрика была читаемой
  }

  void SetUp() override {
    CRSMatrix matrix_a = GenerateBandMatrix(kMatrixSize, kMatrixBandwidth, 2.0);
    CRSMatrix matrix_b = GenerateBandMatrix(kMatrixSize, kMatrixBandwidth, 3.0);
    input_data_ = std::make_tuple(matrix_a, matrix_b);

    expected_output_.rows = kMatrixSize;
    expected_output_.cols = kMatrixSize;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.rows == expected_output_.rows && output_data.cols == expected_output_.cols &&
           !output_data.values.empty();
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  CRSMatrix expected_output_;
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
