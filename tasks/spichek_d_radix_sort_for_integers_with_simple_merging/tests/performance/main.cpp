#include <gtest/gtest.h>

#include <algorithm>
#include <random>

#include "spichek_d_radix_sort_for_integers_with_simple_merging/common/include/common.hpp"
#include "spichek_d_radix_sort_for_integers_with_simple_merging/seq/include/ops_seq.hpp"
#include "task/include/task.hpp"
#include "util/include/perf_test_util.hpp"

namespace spichek_d_radix_sort_for_integers_with_simple_merging {

class SpichekDRadixSortOMP : public SpichekDRadixSortSEQ {
 public:
  using SpichekDRadixSortSEQ::SpichekDRadixSortSEQ;
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
};

class SpichekDRadixSortTBB : public SpichekDRadixSortSEQ {
 public:
  using SpichekDRadixSortSEQ::SpichekDRadixSortSEQ;
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
};

class SpichekDRadixSortSTL : public SpichekDRadixSortSEQ {
 public:
  using SpichekDRadixSortSEQ::SpichekDRadixSortSEQ;
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
};

class SpichekDRadixSortALL : public SpichekDRadixSortSEQ {
 public:
  using SpichekDRadixSortSEQ::SpichekDRadixSortSEQ;
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
};

class SpichekDRadixSortRunPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  const int k_count = 1000000;
  InType input_data;

  void SetUp() override {
    input_data.resize(k_count);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(-100000, 100000);

    for (int i = 0; i < k_count; ++i) {
      input_data[i] = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::ranges::is_sorted(output_data);
  }

  InType GetTestInputData() final {
    return input_data;
  }
};

TEST_P(SpichekDRadixSortRunPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, SpichekDRadixSortALL, SpichekDRadixSortOMP, SpichekDRadixSortSEQ,
                                SpichekDRadixSortSTL, SpichekDRadixSortTBB>(
        PPC_SETTINGS_spichek_d_radix_sort_for_integers_with_simple_merging);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = SpichekDRadixSortRunPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, SpichekDRadixSortRunPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace spichek_d_radix_sort_for_integers_with_simple_merging
