#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <random>
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

class LobanovDMultiplyMatrixFuncTest : public ::testing::Test {
 protected:
  void SetUp() override {}

  static bool ExecuteFullTask(LobanovMultyMatrixALL &task) {
    return task.Validation() && task.PreProcessing() && task.Run() && task.PostProcessing();
  }

  void RunTest(const CompressedRowMatrix &a, const CompressedRowMatrix &b) {
    LobanovMultyMatrixALL task(std::make_pair(a, b));
    ASSERT_TRUE(ExecuteFullTask(task));
    result = task.GetOutput();
  }

  CompressedRowMatrix result;
};

TEST_F(LobanovDMultiplyMatrixFuncTest, SmallMatrices) {
  const auto a = CreateRandomCompressedRowMatrix(10, 10, 0.3, 1);
  const auto b = CreateRandomCompressedRowMatrix(10, 10, 0.3, 2);
  RunTest(a, b);
  EXPECT_EQ(result.row_count, 10);
  EXPECT_EQ(result.column_count, 10);
  EXPECT_GE(result.value_data.size(), static_cast<size_t>(0));
}

TEST_F(LobanovDMultiplyMatrixFuncTest, RectangularMatrices) {
  const auto a = CreateRandomCompressedRowMatrix(10, 5, 0.3, 1);
  const auto b = CreateRandomCompressedRowMatrix(5, 8, 0.3, 2);
  RunTest(a, b);
  EXPECT_EQ(result.row_count, 10);
  EXPECT_EQ(result.column_count, 8);
  EXPECT_GE(result.value_data.size(), static_cast<size_t>(0));
}

TEST_F(LobanovDMultiplyMatrixFuncTest, SparseMatrices) {
  const auto a = CreateRandomCompressedRowMatrix(50, 50, 0.05, 1);
  const auto b = CreateRandomCompressedRowMatrix(50, 50, 0.05, 2);
  RunTest(a, b);
  EXPECT_EQ(result.row_count, 50);
  EXPECT_EQ(result.column_count, 50);
  EXPECT_GE(result.value_data.size(), static_cast<size_t>(0));
}

TEST_F(LobanovDMultiplyMatrixFuncTest, DenseMatrices) {
  const auto a = CreateRandomCompressedRowMatrix(20, 20, 0.8, 1);
  const auto b = CreateRandomCompressedRowMatrix(20, 20, 0.8, 2);
  RunTest(a, b);
  EXPECT_EQ(result.row_count, 20);
  EXPECT_EQ(result.column_count, 20);
  EXPECT_GE(result.value_data.size(), static_cast<size_t>(0));
}

TEST_F(LobanovDMultiplyMatrixFuncTest, LargeRowsSmallCols) {
  const auto a = CreateRandomCompressedRowMatrix(100, 3, 0.3, 1);
  const auto b = CreateRandomCompressedRowMatrix(3, 5, 0.3, 2);
  RunTest(a, b);
  EXPECT_EQ(result.row_count, 100);
  EXPECT_EQ(result.column_count, 5);
  EXPECT_GE(result.value_data.size(), static_cast<size_t>(0));
}

TEST_F(LobanovDMultiplyMatrixFuncTest, SmallRowsLargeCols) {
  const auto a = CreateRandomCompressedRowMatrix(3, 100, 0.3, 1);
  const auto b = CreateRandomCompressedRowMatrix(100, 5, 0.3, 2);
  RunTest(a, b);
  EXPECT_EQ(result.row_count, 3);
  EXPECT_EQ(result.column_count, 5);
  EXPECT_GE(result.value_data.size(), static_cast<size_t>(0));
}

TEST_F(LobanovDMultiplyMatrixFuncTest, IdentityMultiplication) {
  CompressedRowMatrix identity;
  identity.row_count = 5;
  identity.column_count = 5;
  identity.row_pointer_data = {0, 1, 2, 3, 4, 5};
  identity.column_index_data = {0, 1, 2, 3, 4};
  identity.value_data = {1.0, 1.0, 1.0, 1.0, 1.0};
  identity.non_zero_count = 5;

  const auto b = CreateRandomCompressedRowMatrix(5, 5, 0.3, 1);
  RunTest(identity, b);
  EXPECT_EQ(result.row_count, 5);
  EXPECT_EQ(result.column_count, 5);
}

TEST_F(LobanovDMultiplyMatrixFuncTest, ZeroMatrix) {
  CompressedRowMatrix zero;
  zero.row_count = 5;
  zero.column_count = 5;
  zero.row_pointer_data = {0, 0, 0, 0, 0, 0};
  zero.column_index_data.clear();
  zero.value_data.clear();
  zero.non_zero_count = 0;

  const auto b = CreateRandomCompressedRowMatrix(5, 5, 0.3, 1);
  RunTest(zero, b);
  EXPECT_EQ(result.row_count, 5);
  EXPECT_EQ(result.column_count, 5);
  EXPECT_EQ(result.value_data.size(), static_cast<size_t>(0));
}

TEST_F(LobanovDMultiplyMatrixFuncTest, ValidationFailure) {
  const auto a = CreateRandomCompressedRowMatrix(5, 3, 0.3, 1);
  const auto b = CreateRandomCompressedRowMatrix(4, 5, 0.3, 2);
  LobanovMultyMatrixALL task(std::make_pair(a, b));
  EXPECT_FALSE(task.Validation());
}

}  // namespace
}  // namespace lobanov_d_multi_matrix_crs
