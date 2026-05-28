#include <gtest/gtest.h>

#include <cstddef>
#include <random>
#include <tuple>

#include "akhmetov_daniil_strassen_dense_double/all/include/ops_all.hpp"
#include "akhmetov_daniil_strassen_dense_double/common/include/common.hpp"
#include "akhmetov_daniil_strassen_dense_double/omp/include/ops_omp.hpp"
#include "akhmetov_daniil_strassen_dense_double/seq/include/ops_seq.hpp"
#include "akhmetov_daniil_strassen_dense_double/stl/include/ops_stl.hpp"
#include "akhmetov_daniil_strassen_dense_double/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace akhmetov_daniil_strassen_dense_double {

namespace {

class AkhmetovDaniilRunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    TestType n = 1024;

    input_data_.resize(1 + (2 * n * n));
    input_data_.at(0) = static_cast<double>(n);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-10.0, 10.0);

    for (size_t i = 1; i < input_data_.size(); ++i) {
      input_data_.at(i) = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) override {
    if (output_data.empty()) {
      return false;
    }

    size_t n = format::GetN(input_data_);
    return output_data.size() == n * n;
  }

  InType GetTestInputData() override {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(AkhmetovDaniilRunPerfTests, StrassenTestPerf) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks = std::tuple_cat(ppc::util::MakeAllPerfTasks<InType, AkhmetovDStrassenDenseDoubleSEQ>(
                                              PPC_SETTINGS_akhmetov_daniil_strassen_dense_double),
                                          ppc::util::MakeAllPerfTasks<InType, AkhmetovDStrassenDenseDoubleOMP>(
                                              PPC_SETTINGS_akhmetov_daniil_strassen_dense_double),
                                          ppc::util::MakeAllPerfTasks<InType, AkhmetovDStrassenDenseDoubleTBB>(
                                              PPC_SETTINGS_akhmetov_daniil_strassen_dense_double),
                                          ppc::util::MakeAllPerfTasks<InType, AkhmetovDStrassenDenseDoubleSTL>(
                                              PPC_SETTINGS_akhmetov_daniil_strassen_dense_double),
                                          ppc::util::MakeAllPerfTasks<InType, AkhmetovDStrassenDenseDoubleALL>(
                                              PPC_SETTINGS_akhmetov_daniil_strassen_dense_double));

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = AkhmetovDaniilRunPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunStrassenPerfTests, AkhmetovDaniilRunPerfTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace akhmetov_daniil_strassen_dense_double
