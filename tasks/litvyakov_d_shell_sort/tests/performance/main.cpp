#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <random>

#include "litvyakov_d_shell_sort/common/include/common.hpp"
#include "litvyakov_d_shell_sort/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace litvyakov_d_shell_sort {

class LitvyakovDShellSortRunPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    int dim = 1000000;

    InType &in = input_data_;
    in.clear();
    in.reserve(dim);

    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_int_distribution<int> dist(std::numeric_limits<int>::min(), std::numeric_limits<int>::max());

    for (int i = 0; i < dim; ++i) {
      in.push_back(dist(rng));
    }

    test_result_ = in;
    std::ranges::sort(test_result_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != test_result_.size()) {
      return false;
    }
    for (std::size_t i = 0; i < output_data.size(); ++i) {
      if (output_data[i] != test_result_[i]) {
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
  OutType test_result_;
};

TEST_P(LitvyakovDShellSortRunPerfTest, PerfSortTest) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, LitvyakovDShellSortSEQ>(PPC_SETTINGS_litvyakov_d_shell_sort);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = LitvyakovDShellSortRunPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(shellSortPerfTests, LitvyakovDShellSortRunPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace litvyakov_d_shell_sort
