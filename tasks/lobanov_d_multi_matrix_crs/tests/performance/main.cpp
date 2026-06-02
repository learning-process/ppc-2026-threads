#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "lobanov_d_multi_matrix_crs/all/include/ops_all.hpp"
#include "lobanov_d_multi_matrix_crs/common/include/common.hpp"

namespace lobanov_d_multi_matrix_crs {
namespace {

CompressedRowMatrix CreateRandomCompressedRowMatrix(int rows, int cols, double density, int seed = 42) {
  CompressedRowMatrix mat;
  mat.row_count = rows;
  mat.column_count = cols;
  mat.row_pointer_data.clear();
  mat.column_index_data.clear();
  mat.value_data.clear();

  if (rows <= 0 || cols <= 0) {
    mat.row_pointer_data.assign(static_cast<size_t>(rows) + 1, 0);
    return mat;
  }

  density = std::clamp(density, 0.0, 1.0);
  std::mt19937 rng(static_cast<std::mt19937::result_type>(seed));
  std::uniform_real_distribution<double> val_dist(0.1, 10.0);
  std::uniform_real_distribution<double> prob_dist(0.0, 1.0);

  std::vector<std::vector<int>> col_indices_per_row(static_cast<size_t>(rows));
  std::vector<std::vector<double>> values_per_row(static_cast<size_t>(rows));
  int total_nnz = 0;

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (prob_dist(rng) < density) {
        col_indices_per_row[static_cast<size_t>(i)].push_back(j);
        values_per_row[static_cast<size_t>(i)].push_back(val_dist(rng));
        ++total_nnz;
      }
    }
  }

  mat.value_data.reserve(static_cast<size_t>(total_nnz));
  mat.column_index_data.reserve(static_cast<size_t>(total_nnz));
  mat.row_pointer_data.assign(static_cast<size_t>(rows) + 1, 0);

  int offset = 0;
  for (int i = 0; i < rows; ++i) {
    auto &row_cols = col_indices_per_row[static_cast<size_t>(i)];
    auto &row_vals = values_per_row[static_cast<size_t>(i)];

    std::vector<std::pair<int, double>> sorted;
    sorted.reserve(row_cols.size());
    for (size_t k = 0; k < row_cols.size(); ++k) {
      sorted.emplace_back(row_cols[k], row_vals[k]);
    }
    std::ranges::sort(sorted);

    for (const auto &p : sorted) {
      mat.column_index_data.push_back(p.first);
      mat.value_data.push_back(p.second);
    }
    offset += static_cast<int>(row_cols.size());
    mat.row_pointer_data[static_cast<size_t>(i) + 1] = offset;
  }

  if (!mat.row_pointer_data.empty()) {
    mat.row_pointer_data.back() = offset;
  }

  mat.non_zero_count = static_cast<int>(mat.value_data.size());

  return mat;
}

class LobanovDMultiplyMatrixPerfTest : public ::testing::TestWithParam<std::tuple<int, double, std::string>> {
 protected:
  void SetUp() override {
    const auto &params = GetParam();
    dimension_ = std::get<0>(params);
    density_ = std::get<1>(params);
    test_name_ = std::get<2>(params);
  }

  [[nodiscard]] std::pair<CompressedRowMatrix, CompressedRowMatrix> PrepareMatrices() const {
    return {CreateRandomCompressedRowMatrix(dimension_, dimension_, density_, 100),
            CreateRandomCompressedRowMatrix(dimension_, dimension_, density_, 200)};
  }

  static bool ExecuteTask(LobanovMultyMatrixALL &task) {
    return task.Validation() && task.PreProcessing() && task.Run() && task.PostProcessing();
  }

  template <typename Func>
  auto MeasureTime(Func &&func) const {
    const auto start = std::chrono::high_resolution_clock::now();
    std::forward<Func>(func)();
    const auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  }

  void ValidateResult(const CompressedRowMatrix &result) const {
    EXPECT_EQ(result.row_count, dimension_);
    EXPECT_EQ(result.column_count, dimension_);
  }

  void PrintResults(const CompressedRowMatrix &a, const CompressedRowMatrix &b, const CompressedRowMatrix &result,
                    auto duration) const {
    std::cout << "Test: " << test_name_ << '\n';
    std::cout << "Matrix size: " << dimension_ << "x" << dimension_ << '\n';
    std::cout << "Density: " << density_ << '\n';
    std::cout << "NNZ in A: " << a.value_data.size() << '\n';
    std::cout << "NNZ in B: " << b.value_data.size() << '\n';
    std::cout << "NNZ in result: " << result.value_data.size() << '\n';
    std::cout << "Execution time: " << duration.count() << " ms\n";
    std::cout << "------------------------\n";
  }

  void RunPerformanceTest() {
    const auto [a, b] = PrepareMatrices();
    LobanovMultyMatrixALL task(std::make_pair(a, b));
    const auto duration = MeasureTime([&]() { ASSERT_TRUE(ExecuteTask(task)); });
    const auto result = task.GetOutput();
    ValidateResult(result);
    PrintResults(a, b, result, duration);
  }

 private:
  int dimension_ = 0;
  double density_ = 0.0;
  std::string test_name_;
};

TEST_P(LobanovDMultiplyMatrixPerfTest, PerformanceTest) {
  RunPerformanceTest();
}

INSTANTIATE_TEST_SUITE_P(MatrixMultiplicationPerfTests, LobanovDMultiplyMatrixPerfTest,
                         ::testing::Values(std::make_tuple(100, 0.2, "Small_100x100_dense_20%"),
                                           std::make_tuple(500, 0.1, "Medium_500x500_sparse_10%"),
                                           std::make_tuple(1000, 0.05, "Large_1000x1000_sparse_5%"),
                                           std::make_tuple(2000, 0.02, "ExtraLarge_2000x2000_sparse_2%")));

}  // namespace
}  // namespace lobanov_d_multi_matrix_crs
