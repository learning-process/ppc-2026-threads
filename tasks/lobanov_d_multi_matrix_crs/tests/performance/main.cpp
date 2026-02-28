#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <functional>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "lobanov_d_multi_matrix_crs/common/include/common.hpp"
#include "lobanov_d_multi_matrix_crs/seq/include/ops_seq.hpp"

namespace lobanov_d_multi_matrix_crs {

CompressedRowMatrix CreateRandomCompressedRowMatrix(int row_count, int column_count, double density_factor,
                                                    int seed = 42) {
  CompressedRowMatrix result_matrix;
  result_matrix.row_count = row_count;
  result_matrix.column_count = column_count;
  result_matrix.non_zero_count = 0;

  result_matrix.value_data.clear();
  result_matrix.column_index_data.clear();
  result_matrix.row_pointer_data.clear();

  if (row_count <= 0 || column_count <= 0) {
    result_matrix.row_pointer_data.assign(static_cast<std::size_t>(row_count) + 1U, 0);
    return result_matrix;
  }

  density_factor = std::clamp(density_factor, 0.0, 1.0);

  std::mt19937 rng(static_cast<std::mt19937::result_type>(seed));

  std::hash<std::string> hasher;
  std::string param_hash =
      std::to_string(row_count) + "_" + std::to_string(column_count) + "_" + std::to_string(density_factor);
  auto hash_value = static_cast<std::mt19937::result_type>(hasher(param_hash));
  rng.seed(static_cast<std::mt19937::result_type>(seed) + hash_value);

  std::uniform_real_distribution<double> val_dist(0.1, 10.0);
  std::uniform_real_distribution<double> prob_dist(0.0, 1.0);

  std::vector<std::vector<int>> col_indices_per_row(static_cast<std::size_t>(row_count));
  std::vector<std::vector<double>> values_per_row(static_cast<std::size_t>(row_count));

  int nnz_counter = 0;

  for (int i = 0; i < row_count; ++i) {
    for (int j = 0; j < column_count; ++j) {
      if (prob_dist(rng) < density_factor) {
        col_indices_per_row[static_cast<std::size_t>(i)].push_back(j);
        values_per_row[static_cast<std::size_t>(i)].push_back(val_dist(rng));
        ++nnz_counter;
      }
    }
  }

  result_matrix.non_zero_count = nnz_counter;

  if (nnz_counter > 0) {
    result_matrix.value_data.reserve(static_cast<std::size_t>(nnz_counter));
    result_matrix.column_index_data.reserve(static_cast<std::size_t>(nnz_counter));
  }

  result_matrix.row_pointer_data.assign(static_cast<std::size_t>(row_count) + 1U, 0);

  int offset = 0;
  result_matrix.row_pointer_data[0] = 0;

  for (int i = 0; i < row_count; ++i) {
    auto &row_cols = col_indices_per_row[static_cast<std::size_t>(i)];
    auto &row_vals = values_per_row[static_cast<std::size_t>(i)];

    std::vector<std::pair<int, double>> sorted_pairs;
    for (std::size_t k = 0; k < row_cols.size(); ++k) {
      sorted_pairs.emplace_back(row_cols[k], row_vals[k]);
    }
    std::sort(sorted_pairs.begin(), sorted_pairs.end());

    for (const auto &pair : sorted_pairs) {
      result_matrix.column_index_data.push_back(pair.first);
      result_matrix.value_data.push_back(pair.second);
    }

    offset += static_cast<int>(row_cols.size());
    result_matrix.row_pointer_data[static_cast<std::size_t>(i) + 1U] = offset;
  }

  result_matrix.non_zero_count = static_cast<int>(result_matrix.value_data.size());

  if (!result_matrix.row_pointer_data.empty()) {
    result_matrix.row_pointer_data.back() = result_matrix.non_zero_count;
  }

  return result_matrix;
}

class LobanovDMultiplyMatrixPerfTest : public ::testing::TestWithParam<std::tuple<int, double, std::string>> {
 protected:
  void SetUp() override {
    auto params = GetParam();
    dimension_ = std::get<0>(params);
    density_ = std::get<1>(params);
    test_name_ = std::get<2>(params);
  }

  void RunPerformanceTest() {
    auto matrix_a = CreateRandomCompressedRowMatrix(dimension_, dimension_, density_, 100);
    auto matrix_b = CreateRandomCompressedRowMatrix(dimension_, dimension_, density_, 200);

    LobanovMultyMatrixSEQ task(std::make_pair(matrix_a, matrix_b));

    ASSERT_TRUE(task.Validation());

    ASSERT_TRUE(task.PreProcessing());

    auto start = std::chrono::high_resolution_clock::now();

    ASSERT_TRUE(task.Run());

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    ASSERT_TRUE(task.PostProcessing());

    auto result = task.GetOutput();

    EXPECT_EQ(result.row_count, dimension_);
    EXPECT_EQ(result.column_count, dimension_);

    std::cout << "Test: " << test_name_ << std::endl;
    std::cout << "Matrix size: " << dimension_ << "x" << dimension_ << std::endl;
    std::cout << "Density: " << density_ << std::endl;
    std::cout << "NNZ in A: " << matrix_a.non_zero_count << std::endl;
    std::cout << "NNZ in B: " << matrix_b.non_zero_count << std::endl;
    std::cout << "NNZ in result: " << result.non_zero_count << std::endl;
    std::cout << "Execution time: " << duration.count() << " ms" << std::endl;
    std::cout << "------------------------" << std::endl;
  }

  int dimension_;
  double density_;
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
}  // namespace lobanov_d_multi_matrix_crs
