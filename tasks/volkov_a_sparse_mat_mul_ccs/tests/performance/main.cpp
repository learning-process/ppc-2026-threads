#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <tuple>
#include <vector>

#include "util/include/perf_test_util.hpp"
#include "volkov_a_sparse_mat_mul_ccs/common/include/common.hpp"
#include "volkov_a_sparse_mat_mul_ccs/seq/include/ops_seq.hpp"

namespace volkov_a_sparse_mat_mul_ccs {

class VolkovAPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  static SparseMatCCS GenerateRandomSparseMatrix(int rows, int cols) {
    SparseMatCCS mat;
    mat.rows_count = rows;
    mat.cols_count = cols;
    mat.col_ptrs.assign(cols + 1, 0);

    size_t elems_per_col = 50;
    int band = 150;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> val_dist(-5.0, 5.0);

    for (int j = 0; j < cols; ++j) {
      mat.col_ptrs[j] = static_cast<int>(mat.row_indices.size());
      
      int min_row = std::max(0, j - band);
      int max_row = std::min(rows - 1, j + band);
      std::uniform_int_distribution<int> row_dist(min_row, max_row);

      std::vector<int> col_rows;
      col_rows.push_back(j % rows);

      while (col_rows.size() < elems_per_col) {
        int r = row_dist(gen);
        bool exists = false;
        for (int existing_r : col_rows) {
          if (existing_r == r) {
            exists = true;
            break;
          }
        }
        if (!exists) {
          col_rows.push_back(r);
        }
      }

      std::ranges::sort(col_rows);
      
      for (int r : col_rows) {
        mat.row_indices.push_back(r);
        mat.values.push_back(val_dist(gen));
      }
    }
    
    mat.col_ptrs[cols] = static_cast<int>(mat.row_indices.size());
    mat.non_zeros = static_cast<int>(mat.row_indices.size());
    return mat;
  }

  void SetUp() override {
    int n = 20000;
    int k = 20000;
    int m = 20000;
    input_data_ = std::make_tuple(GenerateRandomSparseMatrix(n, k), GenerateRandomSparseMatrix(k, m));
  }

  bool CheckTestOutputData(OutType &matrix_c) final {
    const auto &matrix_a = std::get<0>(input_data_);
    const auto &matrix_b = std::get<1>(input_data_);

    std::vector<double> vector_x(matrix_b.cols_count, 1.0);
    for (int i = 0; i < matrix_b.cols_count; i++) {
      vector_x[i] = (i % 10) * 0.1;
    }

    std::vector<double> vector_bx(matrix_b.rows_count, 0.0);
    for (int j = 0; j < matrix_b.cols_count; ++j) {
      for (int k = matrix_b.col_ptrs[j]; k < matrix_b.col_ptrs[j + 1]; ++k) {
        vector_bx[matrix_b.row_indices[k]] += matrix_b.values[k] * vector_x[j];
      }
    }

    std::vector<double> vector_a_bx(matrix_a.rows_count, 0.0);
    for (int j = 0; j < matrix_a.cols_count; ++j) {
      for (int k = matrix_a.col_ptrs[j]; k < matrix_a.col_ptrs[j + 1]; ++k) {
        vector_a_bx[matrix_a.row_indices[k]] += matrix_a.values[k] * vector_bx[j];
      }
    }

    std::vector<double> vector_cx(matrix_c.rows_count, 0.0);
    for (int j = 0; j < matrix_c.cols_count; ++j) {
      for (int k = matrix_c.col_ptrs[j]; k < matrix_c.col_ptrs[j + 1]; ++k) {
        vector_cx[matrix_c.row_indices[k]] += matrix_c.values[k] * vector_x[j];
      }
    }

    for (int i = 0; i < matrix_c.rows_count; ++i) {
      if (std::abs(vector_a_bx[i] - vector_cx[i]) > 1e-8) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(VolkovAPerfTests, RunPerfTest) {
  ExecuteTest(GetParam());
}

namespace {
const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, VolkovASparseMatMulCcsSeq>(PPC_SETTINGS_volkov_a_sparse_mat_mul_ccs);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = VolkovAPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(VolkovPerfSuite, VolkovAPerfTests, kGtestValues, kPerfTestName);
}  // namespace
}  // namespace volkov_a_sparse_mat_mul_ccs