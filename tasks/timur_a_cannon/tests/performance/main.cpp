#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "timur_a_cannon/common/include/common.hpp"
#include "timur_a_cannon/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"
namespace timur_a_cannon {

class TimurACannonPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 1024;
  InType input_data_;
  OutType res_;

  void SetUp() override {
    int size_block = 64;
    std::vector<std::vector<double>> matrix_a(kCount_, std::vector<double>(kCount_, 2.0));
    std::vector<std::vector<double>> matrix_b(kCount_, std::vector<double>(kCount_, 3.0));

    input_data_ = std::make_tuple(size_block, matrix_a, matrix_b);
    res_ = std::vector<std::vector<double>>(kCount_, std::vector<double>(kCount_, 6144.0));
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if ((res_.size() * res_[0].size()) != (output_data.size() * output_data[0].size())) {
      return false;
    }
    for (size_t i = 0; i < res_.size(); i++) {
      for (size_t j = 0; j < res_[0].size(); j++) {
        if (std::abs(res_[i][j] - output_data[i][j]) > 1e-10) {
          return false;
        }
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(TimurACannonPerfTests, MultiplicationMatrixBlockSchemeCannonPerf) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, TimurACannonMatrixMultiplication>(PPC_SETTINGS_timur_a_cannon);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = TimurACannonPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, TimurACannonPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace timur_a_cannon
