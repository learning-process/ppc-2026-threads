#include <gtest/gtest.h>

#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "util/include/perf_test_util.hpp"
#include "util/include/util.hpp"
#include "vasiliev_m_shell_sort_batcher_merge/common/include/common.hpp"
#include "vasiliev_m_shell_sort_batcher_merge/seq/include/ops_seq.hpp"

namespace vasiliev_m_shell_sort_batcher_merge {

class VasilievMShellSortBatcherMergePerfTestsSEQ : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    std::string abs_path =
        ppc::util::GetAbsoluteTaskPath(PPC_ID_vasiliev_m_shell_sort_batcher_merge, "test_vec_perf.txt");

    std::ifstream file(abs_path);
    if (!file.is_open()) {
      throw std::runtime_error("Wrong path.");
    }

    size_t size = 0;
    file >> size;

    std::vector<int> vec(size);

    for (size_t i = 0; i < size; i++) {
      file >> vec[i];
    }

    input_data_ = vec;
    std::vector<int> basic_vec = input_data_;

    for (int i = 0; i < 500; i++) {
      input_data_.insert(input_data_.end(), basic_vec.begin(), basic_vec.end());
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    (void)output_data;
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(VasilievMShellSortBatcherMergePerfTestsSEQ, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, VasilievMShellSortBatcherMergeSEQ>(
    PPC_SETTINGS_vasiliev_m_shell_sort_batcher_merge);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = VasilievMShellSortBatcherMergePerfTestsSEQ::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerfBatcherTest, VasilievMShellSortBatcherMergePerfTestsSEQ, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace vasiliev_m_shell_sort_batcher_merge
