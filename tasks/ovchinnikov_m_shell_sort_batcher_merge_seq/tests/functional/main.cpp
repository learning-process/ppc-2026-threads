#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>

#include "ovchinnikov_m_shell_sort_batcher_merge_seq/common/include/common.hpp"
#include "ovchinnikov_m_shell_sort_batcher_merge_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace ovchinnikov_m_shell_sort_batcher_merge_seq {

class OvchinnikovMRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    test_data_ = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
  }

  InType GetTestInputData() final {
    return std::get<0>(test_data_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == std::get<1>(test_data_);
  }

 private:
  TestType test_data_;
};

namespace {

TEST_P(OvchinnikovMRunFuncTestsThreads, ShellSortBatcherMergeFuncTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 8> kTestParam = {
    TestType{{}, {}, "Empty"},
    TestType{{42}, {42}, "Single"},
    TestType{{2, 1}, {1, 2}, "TwoElements"},
    TestType{{5, 4, 3, 2, 1}, {1, 2, 3, 4, 5}, "Reverse"},
    TestType{{1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}, "Sorted"},
    TestType{{0, -1, 5, -10, 3}, {-10, -1, 0, 3, 5}, "MixedSigns"},
    TestType{{7, 7, 7, 7}, {7, 7, 7, 7}, "Duplicates"},
    TestType{{9, 0, 8, 1, 7, 2, 6, 3}, {0, 1, 2, 3, 6, 7, 8, 9}, "EvenOddLength"}};

const auto kTestTasksList = ppc::util::AddFuncTask<OvchinnikovMShellSortBatcherMergeSEQ, InType>(
    kTestParam, PPC_SETTINGS_ovchinnikov_m_shell_sort_batcher_merge_seq);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = OvchinnikovMRunFuncTestsThreads::PrintFuncTestName<OvchinnikovMRunFuncTestsThreads>;

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
INSTANTIATE_TEST_SUITE_P(ShellSortBatcherMergeTests, OvchinnikovMRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace ovchinnikov_m_shell_sort_batcher_merge_seq
