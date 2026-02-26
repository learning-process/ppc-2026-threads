#include <gtest/gtest.h>

#include <cstddef>
#include <random>
#include <set>
#include <tuple>

#include "maslova_u_mult_matr_crs/common/include/common.hpp"
#include "maslova_u_mult_matr_crs/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace maslova_u_mult_matr_crs {

namespace {

CRSMatrix CreateUniqueRandomCRS(int rows, int cols, double density) {
  CRSMatrix matrix;
  matrix.rows = rows;
  matrix.cols = cols;
  matrix.row_ptr.assign(rows + 1, 0);

  std::mt19937 generator(std::random_device{}());
  std::uniform_real_distribution<double> val_dist(-100.0, 100.0);
  std::uniform_int_distribution<int> col_dist(0, cols - 1);

  int elements_per_row = static_cast<int>(cols * density);
  if (elements_per_row < 1 && density > 0) {
    elements_per_row = 1;
  }

  for (int i = 0; i < rows; ++i) {
    std::set<int> selected_columns;
    while (selected_columns.size() < static_cast<size_t>(elements_per_row)) {
      selected_columns.insert(col_dist(generator));
    }

    for (int col : selected_columns) {
      matrix.col_ind.push_back(col);
      matrix.values.push_back(val_dist(generator));
    }
    matrix.row_ptr[i + 1] = static_cast<int>(matrix.values.size());
  }

  return matrix;
}

}  // namespace

class MaslovaUMultMatrRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    int size = 1000;
    double density = 0.01;
    CRSMatrix a = CreateUniqueRandomCRS(size, size, density);
    CRSMatrix b = CreateUniqueRandomCRS(size, size, density);
    input_data_ = std::make_tuple(a, b);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const auto &matrix_a = std::get<0>(input_data_);
    const auto &matrix_b = std::get<1>(input_data_);

    return std::tie(output_data.rows, output_data.cols) == std::tie(matrix_a.rows, matrix_b.cols);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(MaslovaUMultMatrRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, MaslovaUMultMatrSEQ>(PPC_SETTINGS_maslova_u_mult_matr_crs);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = MaslovaUMultMatrRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, MaslovaUMultMatrRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace maslova_u_mult_matr_crs
