#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <tuple>
#include <vector>

#include "kulik_a_mat_mul_double_ccs/common/include/common.hpp"
#include "kulik_a_mat_mul_double_ccs/omp/include/ops_omp.hpp"
#include "kulik_a_mat_mul_double_ccs/seq/include/ops_seq.hpp"
#include "kulik_a_mat_mul_double_ccs/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace kulik_a_mat_mul_double_ccs {

class KulikARunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    size_t size = 20000;
    CCS &a = std::get<0>(input_data_);
    CCS &b = std::get<1>(input_data_);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist_val(-10.0, 10.0);
    a.m = size;
    a.n = size;
    b.m = size;
    b.n = size;
    a.col_ind.assign(a.n + 1, 0);
    b.col_ind.assign(b.n + 1, 0);
    for (size_t j = 0; j < a.n; ++j) {
      size_t left = (j > 50) ? (j - 50) : 0;
      size_t right = std::min(a.n, j + 50);
      a.col_ind[j + 1] = a.col_ind[j] + right - left;
      for (size_t k = left; k < right; ++k) {
        double r = dist_val(gen);
        a.row.push_back(k);
        a.value.push_back(r);
      }
    }
    a.nz = a.row.size();
    for (size_t j = 0; j < b.n; ++j) {
      size_t left = (j > 50) ? (j - 50) : 0;
      size_t right = std::min(b.n, j + 50);
      b.col_ind[j + 1] = b.col_ind[j] + right - left;
      for (size_t k = left; k < right; ++k) {
        double r = dist_val(gen);
        b.row.push_back(k);
        b.value.push_back(r);
      }
    }
    b.nz = b.row.size();

    input_data_ = std::make_tuple(a, b);
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
    ppc::util::MakeAllPerfTasks<InType, KulikAMatMulDoubleCcsSEQ, KulikAMatMulDoubleCcsOMP, KulikAMatMulDoubleCcsTBB>(
        PPC_SETTINGS_kulik_a_mat_mul_double_ccs);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KulikARunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KulikARunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kulik_a_mat_mul_double_ccs
