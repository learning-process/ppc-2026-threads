#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "muhammadkhon_i_stressen_alg/common/include/common.hpp"
#include "muhammadkhon_i_stressen_alg/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace muhammadkhon_i_stressen_alg {

class MuhammadkhonIStressenAlgPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  void SetUp() override {
    const size_t rc = 512;
    const size_t size = rc * rc;

    input_data_.a_rows = rc;
    input_data_.a_cols_b_rows = rc;
    input_data_.b_cols = rc;

    input_data_.a.assign(size, 0.0);
    input_data_.b.assign(size, 0.0);

    for (size_t i = 0; i < size; i++) {
      input_data_.a[i] = static_cast<double>(i % 100);
      input_data_.b[i] = static_cast<double>((i + 1) % 100);
    }

    expected_output_.assign(size, 0.0);

    for (size_t i = 0; i < rc; ++i) {
      for (size_t k = 0; k < rc; ++k) {
        double temp = input_data_.a[(i * rc) + k];
        for (size_t j = 0; j < rc; ++j) {
          expected_output_[(i * rc) + j] += temp * input_data_.b[(k * rc) + j];
        }
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (expected_output_.size() != output_data.size()) {
      return false;
    }
    constexpr double kEpsilon = 1e-9;
    for (size_t i = 0; i < expected_output_.size(); ++i) {
      if (std::abs(expected_output_[i] - output_data[i]) > kEpsilon) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  std::vector<double> expected_output_;
};

TEST_P(MuhammadkhonIStressenAlgPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, MuhammadkhonIStressenAlgSEQ>(PPC_SETTINGS_muhammadkhon_i_stressen_alg);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = MuhammadkhonIStressenAlgPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, MuhammadkhonIStressenAlgPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace muhammadkhon_i_stressen_alg
