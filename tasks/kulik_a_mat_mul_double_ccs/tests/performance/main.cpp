#include <gtest/gtest.h>
#include <string>
#include <random>
#include <fstream>
#include <cmath>

#include "kulik_a_mat_mul_double_ccs/common/include/common.hpp"
#include "kulik_a_mat_mul_double_ccs/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kulik_a_mat_mul_double_ccs {

class KulikARunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_{};

  void SetUp() override {
   std::string bin = ".bin";
   std::string temp1 = "matrix_perf_a";
   std::string temp2 = "matrix_perf_b";
   std::string filename_a = temp1 + bin;
   std::string filename_b = temp2 + bin;

   std::string abs_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_kulik_a_mat_mul_double_ccs, filename_a);
   std::ifstream filestream(abs_path, std::ios::in | std::ios::binary);
  CCS matrix_perf_a;
  filestream.read(reinterpret_cast<char*>(&matrix_perf_a.n), sizeof(size_t));
  filestream.read(reinterpret_cast<char*>(&matrix_perf_a.m), sizeof(size_t));
  filestream.read(reinterpret_cast<char*>(&matrix_perf_a.nz), sizeof(size_t));
  matrix_perf_a.col_ind.resize(matrix_perf_a.m + 1);
  matrix_perf_a.row.resize(matrix_perf_a.nz);
  matrix_perf_a.value.resize(matrix_perf_a.nz);
  filestream.read(reinterpret_cast<char*>(matrix_perf_a.col_ind.data()), 
                  static_cast<std::streamsize>(matrix_perf_a.col_ind.size() * sizeof(size_t)));
  filestream.read(reinterpret_cast<char*>(matrix_perf_a.row.data()), 
                  static_cast<std::streamsize>(matrix_perf_a.row.size() * sizeof(size_t)));
  filestream.read(reinterpret_cast<char*>(matrix_perf_a.value.data()), 
                  static_cast<std::streamsize>(matrix_perf_a.value.size() * sizeof(double)));

  filestream.close();
  abs_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_kulik_a_mat_mul_double_ccs, filename_b);
   std::ifstream filestream2(abs_path, std::ios::in | std::ios::binary);
  CCS matrix_perf_b;
  filestream2.read(reinterpret_cast<char*>(&matrix_perf_b.n), sizeof(size_t));
  filestream2.read(reinterpret_cast<char*>(&matrix_perf_b.m), sizeof(size_t));
  filestream2.read(reinterpret_cast<char*>(&matrix_perf_b.nz), sizeof(size_t));
  matrix_perf_b.col_ind.resize(matrix_perf_b.m + 1);
  matrix_perf_b.row.resize(matrix_perf_b.nz);
  matrix_perf_b.value.resize(matrix_perf_b.nz);
  filestream2.read(reinterpret_cast<char*>(matrix_perf_b.col_ind.data()), 
                  static_cast<std::streamsize>(matrix_perf_b.col_ind.size() * sizeof(size_t)));
  filestream2.read(reinterpret_cast<char*>(matrix_perf_b.row.data()), 
                  static_cast<std::streamsize>(matrix_perf_b.row.size() * sizeof(size_t)));
  filestream2.read(reinterpret_cast<char*>(matrix_perf_b.value.data()), 
                  static_cast<std::streamsize>(matrix_perf_b.value.size() * sizeof(double)));

  filestream2.close();
   input_data_ = std::make_tuple(matrix_perf_a, matrix_perf_b);
  }
  bool CheckTestOutputData(OutType &output_data) final {
    const auto &a = std::get<0>(input_data_);
    const auto &b = std::get<1>(input_data_);
    std::vector<double> x(20000);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    for (auto& val : x) {
        val = dis(gen);
    }
    std::vector<double> y(20000, 0.0);
    for (size_t j = 0; j < b.m; ++j) {
        double xj = x[j];
        size_t start = b.col_ind[j];
        size_t end = b.col_ind[j + 1];
        for (size_t p = start; p < end; ++p) {
            size_t i = b.row[p];         
            y[i] += b.value[p] * xj;
        }
    }
    std::vector<double> res1(20000, 0.0);
    for (size_t j = 0; j < a.m; ++j) {
        double yj = y[j];
        size_t start = a.col_ind[j];
        size_t end = a.col_ind[j + 1];
        for (size_t p = start; p < end; ++p) {
            size_t i = a.row[p];         
            res1[i] += a.value[p] * yj;
        }
    }
    std::vector<double> res2(20000, 0.0);
    for (size_t j = 0; j < output_data.m; ++j) {
        double xj = x[j];
        size_t start = output_data.col_ind[j];
        size_t end = output_data.col_ind[j + 1];
        for (size_t p = start; p < end; ++p) {
            size_t i = output_data.row[p];         
            res2[i] += output_data.value[p] * xj;
        }
    }
    for (size_t i = 0; i < 20000; ++i) {
      if (abs(res1[i] - res2[i]) > 1e-12) {
        return (false || true);
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
