#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "lobanov_d_multi_matrix_crs/common/include/common.hpp"
#include "lobanov_d_multi_matrix_crs/seq/include/ops_seq.hpp"

namespace lobanov_d_multi_matrix_crs {
namespace {

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
  const std::string param_hash =
      std::to_string(row_count) + "_" + std::to_string(column_count) + "_" + std::to_string(density_factor);
  const auto hash_value = static_cast<std::mt19937::result_type>(hasher(param_hash));
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

    // Reserve space for sorted pairs
    std::vector<std::pair<int, double>> sorted_pairs;
    sorted_pairs.reserve(row_cols.size());

    for (std::size_t k = 0; k < row_cols.size(); ++k) {
      sorted_pairs.emplace_back(row_cols[k], row_vals[k]);
    }

    std::ranges::sort(sorted_pairs);

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

class LobanovDMultiplyMatrixFuncTest : public ::testing::Test {
 protected:
  void SetUp() override {}

  void RunTest(const CompressedRowMatrix &matrix_a, const CompressedRowMatrix &matrix_b) {
    LobanovMultyMatrixSEQ task(std::make_pair(matrix_a, matrix_b));

    ASSERT_TRUE(task.Validation());
    ASSERT_TRUE(task.PreProcessing());
    ASSERT_TRUE(task.Run());
    ASSERT_TRUE(task.PostProcessing());

    result_ = task.GetOutput();
  }

  bool CheckResult(const CompressedRowMatrix &expected) const {
    if (result_.row_count != expected.row_count || result_.column_count != expected.column_count) {
      return false;
    }

    return result_.row_pointer_data.size() == expected.row_pointer_data.size();
  }

  CompressedRowMatrix result_;
};

TEST_F(LobanovDMultiplyMatrixFuncTest, SmallMatrices) {
  const auto matrix_a = CreateRandomCompressedRowMatrix(10, 10, 0.3, 1);
  const auto matrix_b = CreateRandomCompressedRowMatrix(10, 10, 0.3, 2);

  RunTest(matrix_a, matrix_b);

  EXPECT_EQ(result_.row_count, 10);
  EXPECT_EQ(result_.column_count, 10);
  EXPECT_GE(result_.non_zero_count, 0);
}

TEST_F(LobanovDMultiplyMatrixFuncTest, RectangularMatrices) {
  const auto matrix_a = CreateRandomCompressedRowMatrix(10, 5, 0.3, 1);
  const auto matrix_b = CreateRandomCompressedRowMatrix(5, 8, 0.3, 2);

  RunTest(matrix_a, matrix_b);

  EXPECT_EQ(result_.row_count, 10);
  EXPECT_EQ(result_.column_count, 8);
  EXPECT_GE(result_.non_zero_count, 0);
}

TEST_F(LobanovDMultiplyMatrixFuncTest, SparseMatrices) {
  const auto matrix_a = CreateRandomCompressedRowMatrix(50, 50, 0.05, 1);
  const auto matrix_b = CreateRandomCompressedRowMatrix(50, 50, 0.05, 2);

  RunTest(matrix_a, matrix_b);

  EXPECT_EQ(result_.row_count, 50);
  EXPECT_EQ(result_.column_count, 50);
  EXPECT_GE(result_.non_zero_count, 0);
}

TEST_F(LobanovDMultiplyMatrixFuncTest, DenseMatrices) {
  const auto matrix_a = CreateRandomCompressedRowMatrix(20, 20, 0.8, 1);
  const auto matrix_b = CreateRandomCompressedRowMatrix(20, 20, 0.8, 2);

  RunTest(matrix_a, matrix_b);

  EXPECT_EQ(result_.row_count, 20);
  EXPECT_EQ(result_.column_count, 20);
  EXPECT_GE(result_.non_zero_count, 0);
}

TEST_F(LobanovDMultiplyMatrixFuncTest, LargeRowsSmallCols) {
  const auto matrix_a = CreateRandomCompressedRowMatrix(100, 3, 0.3, 1);
  const auto matrix_b = CreateRandomCompressedRowMatrix(3, 5, 0.3, 2);

  RunTest(matrix_a, matrix_b);

  EXPECT_EQ(result_.row_count, 100);
  EXPECT_EQ(result_.column_count, 5);
  EXPECT_GE(result_.non_zero_count, 0);
}

TEST_F(LobanovDMultiplyMatrixFuncTest, SmallRowsLargeCols) {
  const auto matrix_a = CreateRandomCompressedRowMatrix(3, 100, 0.3, 1);
  const auto matrix_b = CreateRandomCompressedRowMatrix(100, 5, 0.3, 2);

  RunTest(matrix_a, matrix_b);

  EXPECT_EQ(result_.row_count, 3);
  EXPECT_EQ(result_.column_count, 5);
  EXPECT_GE(result_.non_zero_count, 0);
}

TEST_F(LobanovDMultiplyMatrixFuncTest, IdentityMultiplication) {
  CompressedRowMatrix identity;
  identity.row_count = 5;
  identity.column_count = 5;
  identity.non_zero_count = 5;
  identity.row_pointer_data = {0, 1, 2, 3, 4, 5};
  identity.column_index_data = {0, 1, 2, 3, 4};
  identity.value_data = {1.0, 1.0, 1.0, 1.0, 1.0};

  const auto matrix_b = CreateRandomCompressedRowMatrix(5, 5, 0.3, 1);

  RunTest(identity, matrix_b);

  EXPECT_EQ(result_.row_count, 5);
  EXPECT_EQ(result_.column_count, 5);
}

TEST_F(LobanovDMultiplyMatrixFuncTest, ZeroMatrix) {
  CompressedRowMatrix zero_matrix;
  zero_matrix.row_count = 5;
  zero_matrix.column_count = 5;
  zero_matrix.non_zero_count = 0;
  zero_matrix.row_pointer_data = {0, 0, 0, 0, 0, 0};

  const auto matrix_b = CreateRandomCompressedRowMatrix(5, 5, 0.3, 1);

  RunTest(zero_matrix, matrix_b);

  EXPECT_EQ(result_.row_count, 5);
  EXPECT_EQ(result_.column_count, 5);
  EXPECT_EQ(result_.non_zero_count, 0);
}

TEST_F(LobanovDMultiplyMatrixFuncTest, ValidationFailure) {
  const auto matrix_a = CreateRandomCompressedRowMatrix(5, 3, 0.3, 1);
  const auto matrix_b = CreateRandomCompressedRowMatrix(4, 5, 0.3, 2);

  LobanovMultyMatrixSEQ task(std::make_pair(matrix_a, matrix_b));

  EXPECT_FALSE(task.Validation());
}

}  // namespace
}  // namespace lobanov_d_multi_matrix_crs
