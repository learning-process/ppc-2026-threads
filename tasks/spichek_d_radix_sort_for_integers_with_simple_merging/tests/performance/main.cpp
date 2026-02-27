#include <gtest/gtest.h>

#include <algorithm>
#include <random>

#include "spichek_d_radix_sort_for_integers_with_simple_merging/common/include/common.hpp"
#include "spichek_d_radix_sort_for_integers_with_simple_merging/seq/include/ops_seq.hpp"
#include "task/include/task.hpp"  // Добавлено для ppc::task::TypeOfTask
#include "util/include/perf_test_util.hpp"

namespace spichek_d_radix_sort_for_integers_with_simple_merging {

class RadixSortOMP : public RadixSortSEQ {
 public:
  using RadixSortSEQ::RadixSortSEQ;
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
};

class RadixSortTBB : public RadixSortSEQ {
 public:
  using RadixSortSEQ::RadixSortSEQ;
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
};

class RadixSortSTL : public RadixSortSEQ {
 public:
  using RadixSortSEQ::RadixSortSEQ;
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
};

class RadixSortALL : public RadixSortSEQ {
 public:
  using RadixSortSEQ::RadixSortSEQ;
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
};

class RadixSortRunPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  // Исправлено именование согласно readability-identifier-naming
  const int kCount = 1000000;
  InType input_data;  // Удалена избыточная инициализация {}

  void SetUp() override {
    input_data.resize(kCount);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(-100000, 100000);

    for (int i = 0; i < kCount; ++i) {
      input_data[i] = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::ranges::is_sorted(output_data);  // Переход на ranges
  }

  InType GetTestInputData() final {
    return input_data;
  }
};

TEST_P(RadixSortRunPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, RadixSortALL, RadixSortOMP, RadixSortSEQ, RadixSortSTL, RadixSortTBB>(
        PPC_SETTINGS_spichek_d_radix_sort_for_integers_with_simple_merging);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = RadixSortRunPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, RadixSortRunPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace spichek_d_radix_sort_for_integers_with_simple_merging
