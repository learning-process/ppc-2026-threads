#include <gtest/gtest.h>
#include <random>

#include "chetverikova_e_shell_sort_simple_merge/common/include/common.hpp"
#include "chetverikova_e_shell_sort_simple_merge/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace chetverikova_e_shell_sort_simple_merge {

class ChetverikovaERunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;
  OutType expected_data_;

  void SetUp() override {
    constexpr std::size_t kSize = 40;
    std::random_device random_device;
    std::mt19937 generator(random_device());
    std::uniform_int_distribution<int> distribution(0, 100000);

    input_data_.resize(kSize);
    for (auto &value : input_data_) {
      value = distribution(generator);
    }

    expected_data_ = input_data_;
    std::ranges::sort(expected_data_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != input_data_.size()) {
      return false;
    }
    std::cout << "Expected size: " << expected_data_.size() << std::endl;
    std::cout << "Output size: " << output_data.size() << std::endl;
    
    for (size_t i = 0; i < output_data.size(); i++) {
      if (output_data[i] != expected_data_[i]) {
        std::cout << "Mismatch detected!\n";
        std::cout << "Mismatch at index " << i
                    << ": exp=" << expected_data_[i]
                    << ", got=" << output_data[i]
                    << std::endl;
        return false;
      }
    }
    return true;
    //return expected_data_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(ChetverikovaERunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ChetverikovaEShellSortSimpleMergeSEQ>(
    PPC_SETTINGS_chetverikova_e_shell_sort_simple_merge);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ChetverikovaERunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ChetverikovaERunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace chetverikova_e_shell_sort_simple_merge
