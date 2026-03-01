#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <tuple>
#include <vector>

#include "kulik_a_mat_mul_double_ccs/common/include/common.hpp"
#include "kulik_a_mat_mul_double_ccs/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kulik_a_mat_mul_double_ccs {

class KulikARunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    size_t n = 20000;
    size_t m = 20000;
    size_t k = 20000;         // Размеры матриц (n x k) и (k x m)
    size_t nnz_per_col = 50;  // кол-во ненулевых элементов в столбце
    size_t band_width = 175;  // ленточная матрица для лучшей локальности данных

    auto generate_ccs = [&](size_t rows, size_t cols) {
      CCS mat;
      mat.n = rows;
      mat.m = cols;
      mat.col_ind.assign(cols + 1, 0);
      mat.row.reserve(cols * nnz_per_col);
      mat.value.reserve(cols * nnz_per_col);

      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<double> dist_val(-10.0, 10.0);

      for (size_t j = 0; j < cols; ++j) {
        mat.col_ind[j] = mat.row.size();
        size_t min_r = 0;
        if (j >= band_width) {
          min_r = j - band_width;
        }
        size_t max_r = std::min(rows - 1, j + band_width);

        std::vector<size_t> current_rows;
        current_rows.push_back(j % rows);

        std::uniform_int_distribution<size_t> dist_row(min_r, max_r);
        while (current_rows.size() < nnz_per_col) {
          size_t r = dist_row(gen);
          if (std::ranges::find(current_rows, r) == current_rows.end()) {
            current_rows.push_back(r);
          }
        }
        std::ranges::sort(current_rows);
        for (size_t row : current_rows) {
          mat.row.push_back(row);
          mat.value.push_back(dist_val(gen));
        }
      }
      mat.col_ind[cols] = mat.row.size();
      mat.nz = mat.row.size();
      return mat;
    };

    input_data_ = std::make_tuple(generate_ccs(n, k), generate_ccs(k, m));
  }
  bool CheckTestOutputData(OutType &output_data) final {
    const auto &a = std::get<0>(input_data_);
    const auto &b = std::get<1>(input_data_);
    std::vector<double> x(b.m);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    for (auto &val : x) {
      val = dis(gen);
    }
    std::vector<double> y(b.n, 0.0);
    for (size_t j = 0; j < b.m; ++j) {
      double xj = x[j];
      size_t start = b.col_ind[j];
      size_t end = b.col_ind[j + 1];
      for (size_t pc = start; pc < end; ++pc) {
        size_t i = b.row[pc];
        y[i] += b.value[pc] * xj;
      }
    }
    std::vector<double> res1(a.m, 0.0);
    for (size_t j = 0; j < a.m; ++j) {
      double yj = y[j];
      size_t start = a.col_ind[j];
      size_t end = a.col_ind[j + 1];
      for (size_t pc = start; pc < end; ++pc) {
        size_t i = a.row[pc];
        res1[i] += a.value[pc] * yj;
      }
    }
    std::vector<double> res2(output_data.m, 0.0);
    for (size_t j = 0; j < output_data.m; ++j) {
      double xj = x[j];
      size_t start = output_data.col_ind[j];
      size_t end = output_data.col_ind[j + 1];
      for (size_t pc = start; pc < end; ++pc) {
        size_t i = output_data.row[pc];
        res2[i] += output_data.value[pc] * xj;
      }
    }
    for (size_t i = 0; i < output_data.n; ++i) {
      if (std::abs(res1[i] - res2[i]) > 1e-10) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KulikARunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KulikAMatMulDoubleCcsSEQ>(PPC_SETTINGS_kulik_a_mat_mul_double_ccs);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KulikARunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KulikARunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kulik_a_mat_mul_double_ccs
