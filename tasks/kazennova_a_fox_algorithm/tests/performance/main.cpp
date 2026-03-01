#include <gtest/gtest.h>

#include <cstddef>
#include <random>
#include <vector>

#include "kazennova_a_fox_algorithm/common/include/common.hpp"
#include "kazennova_a_fox_algorithm/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kazennova_a_fox_algorithm {

class KazennovaAPerfTestSeq : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kMatrixSize = 500;
  InType input_data_;

  void SetUp() override {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-10.0, 10.0);

    // Матрица A (kMatrixSize × kMatrixSize)
    input_data_.A.rows = kMatrixSize;
    input_data_.A.cols = kMatrixSize;
    input_data_.A.data.resize(static_cast<size_t>(kMatrixSize) * kMatrixSize);
    for (int i = 0; i < kMatrixSize * kMatrixSize; ++i) {
      input_data_.A.data[i] = dis(gen);
    }

    // Матрица B (kMatrixSize × kMatrixSize)
    input_data_.B.rows = kMatrixSize;
    input_data_.B.cols = kMatrixSize;
    input_data_.B.data.resize(static_cast<size_t>(kMatrixSize) * kMatrixSize);
    for (int i = 0; i < kMatrixSize * kMatrixSize; ++i) {
      input_data_.B.data[i] = dis(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.rows == input_data_.A.rows && output_data.cols == input_data_.B.cols &&
           output_data.data.size() == static_cast<size_t>(output_data.rows) * output_data.cols;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KazennovaAPerfTestSeq, RunPerfTests) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KazennovaATestTaskSEQ>(PPC_SETTINGS_kazennova_a_fox_algorithm);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KazennovaAPerfTestSeq::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunPerfTests, KazennovaAPerfTestSeq, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kazennova_a_fox_algorithm
