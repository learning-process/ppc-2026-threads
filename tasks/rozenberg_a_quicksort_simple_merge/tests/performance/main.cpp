#include <gtest/gtest.h>

#include <fstream>
#include <string>

#include "rozenberg_a_quicksort_simple_merge/common/include/common.hpp"
#include "rozenberg_a_quicksort_simple_merge/omp/include/ops_omp.hpp"
#include "rozenberg_a_quicksort_simple_merge/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"
#include "util/include/util.hpp"

namespace rozenberg_a_quicksort_simple_merge {

class RozenbergARunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  void SetUp() override {
    input_data_.clear();
    output_data_.clear();
    std::string abs_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_rozenberg_a_quicksort_simple_merge, "perf_test.txt");
    std::ifstream file(abs_path);

    if (file.is_open()) {
      int size = 0;
      file >> size;

      InType input_data(size);
      for (int i = 0; i < size; i++) {
        file >> input_data[i];
      }

      OutType output_data(size);
      for (int i = 0; i < size; i++) {
        file >> output_data[i];
      }
      input_data_ = input_data;
      output_data_ = output_data;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (output_data_ == output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType output_data_;
};

TEST_P(RozenbergARunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, RozenbergAQuicksortSimpleMergeSEQ, RozenbergAQuicksortSimpleMergeOMP>(
        PPC_SETTINGS_rozenberg_a_quicksort_simple_merge);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = RozenbergARunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, RozenbergARunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace rozenberg_a_quicksort_simple_merge
