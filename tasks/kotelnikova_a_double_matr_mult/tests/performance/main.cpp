#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <random>
#include <utility>
#include <vector>

#include "kotelnikova_a_double_matr_mult/all/include/ops_all.hpp"
#include "kotelnikova_a_double_matr_mult/common/include/common.hpp"
#include "kotelnikova_a_double_matr_mult/omp/include/ops_omp.hpp"
#include "kotelnikova_a_double_matr_mult/seq/include/ops_seq.hpp"
#include "kotelnikova_a_double_matr_mult/stl/include/ops_stl.hpp"
#include "kotelnikova_a_double_matr_mult/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace kotelnikova_a_double_matr_mult {

namespace {
SparseMatrixCCS CreateTestMatrix(int size, double density) {
  SparseMatrixCCS matrix(size, size);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> value_dist(-10.0, 10.0);
  std::uniform_real_distribution<double> density_dist(0.0, 1.0);
  std::vector<std::vector<std::pair<int, double>>> columns(size);

  for (int j = 0; j < size; ++j) {
    for (int i = 0; i < size; ++i) {
      if (density_dist(gen) < density) {
        columns[j].emplace_back(i, value_dist(gen));
      }
    }
    std::sort(columns[j].begin(), columns[j].end());
  }

  matrix.col_ptrs[0] = 0;
  for (int j = 0; j < size; ++j) {
    matrix.col_ptrs[j + 1] = matrix.col_ptrs[j] + static_cast<int>(columns[j].size());
    for (const auto &[row, val] : columns[j]) {
      matrix.row_indices.push_back(row);
      matrix.values.push_back(val);
    }
  }

  return matrix;
}
}  // namespace

class KotelnikovaARunPerfTestSEQ : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kMatrixSize = 700;
  static constexpr double kDensity = 0.1;
  InType input_data_;

  void SetUp() override {
    const SparseMatrixCCS a = CreateTestMatrix(kMatrixSize, kDensity);
    const SparseMatrixCCS b = CreateTestMatrix(kMatrixSize, kDensity);
    input_data_ = std::make_pair(a, b);
  }

  bool CheckTestOutputData(OutType &output_data) override {
    return (output_data.rows == kMatrixSize && output_data.cols == kMatrixSize);
  }

  InType GetTestInputData() override {
    return input_data_;
  }
};

TEST_P(KotelnikovaARunPerfTestSEQ, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KotelnikovaATaskSEQ, KotelnikovaATaskOMP, KotelnikovaATaskTBB,
                                KotelnikovaATaskSTL, KotelnikovaATaskALL>(PPC_SETTINGS_kotelnikova_a_double_matr_mult);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KotelnikovaARunPerfTestSEQ::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(SparseMatrixMultPerfTests, KotelnikovaARunPerfTestSEQ, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kotelnikova_a_double_matr_mult
