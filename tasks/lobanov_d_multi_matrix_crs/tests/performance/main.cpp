#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <tuple>
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

class LobanovDMultiplyMatrixPerfTest : public ::testing::TestWithParam<std::tuple<int, double, std::string>> {
 protected:
  void SetUp() override {
    const auto &params = GetParam();
    dimension = std::get<0>(params);
    density = std::get<1>(params);
    test_name = std::get<2>(params);
  }

  // Функция для подготовки данных
  std::pair<CompressedRowMatrix, CompressedRowMatrix> PrepareMatrices() const {
    return {CreateRandomCompressedRowMatrix(dimension, dimension, density, 100),
            CreateRandomCompressedRowMatrix(dimension, dimension, density, 200)};
  }

  // Функция для выполнения задачи
  bool ExecuteTask(LobanovMultyMatrixSEQ &task) const {
    return task.Validation() && task.PreProcessing() && task.Run() && task.PostProcessing();
  }

  // Функция для измерения времени выполнения
  template <typename Func>
  auto MeasureTime(Func &&func) const {
    const auto start = std::chrono::high_resolution_clock::now();
    func();
    const auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  }

  // Функция для проверки результатов
  void ValidateResult(const CompressedRowMatrix &result) const {
    EXPECT_EQ(result.row_count, dimension);
    EXPECT_EQ(result.column_count, dimension);
  }

  // Функция для вывода результатов
  void PrintResults(const std::string &test_name, int dimension, double density, const CompressedRowMatrix &matrix_a,
                    const CompressedRowMatrix &matrix_b, const CompressedRowMatrix &result, auto duration) const {
    std::cout << "Test: " << test_name << '\n';
    std::cout << "Matrix size: " << dimension << "x" << dimension << '\n';
    std::cout << "Density: " << density << '\n';
    std::cout << "NNZ in A: " << matrix_a.non_zero_count << '\n';
    std::cout << "NNZ in B: " << matrix_b.non_zero_count << '\n';
    std::cout << "NNZ in result: " << result.non_zero_count << '\n';
    std::cout << "Execution time: " << duration.count() << " ms\n";
    std::cout << "------------------------\n";
  }

  // Основная функция теперь очень простая
  void RunPerformanceTest() {
    // Подготовка
    const auto [matrix_a, matrix_b] = PrepareMatrices();
    LobanovMultyMatrixSEQ task(std::make_pair(matrix_a, matrix_b));

    // Выполнение с измерением времени
    const auto duration = MeasureTime([&]() { ASSERT_TRUE(ExecuteTask(task)); });

    // Получение и проверка результата
    const auto result = task.GetOutput();
    ValidateResult(result);

    // Вывод результатов
    PrintResults(test_name, dimension, density, matrix_a, matrix_b, result, duration);
  }

  int dimension = 0;
  double density = 0.0;
  std::string test_name;
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
