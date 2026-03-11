#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>

#include "frolova_s_radix_sort_double/common/include/common.hpp"
#include "frolova_s_radix_sort_double/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace frolova_s_radix_sort_double {

class FrolovaSRadixSortDoubleRunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  InType GetTestInputData() final {
    return input_data;
  }

 protected:
  InType input_data;
  OutType expected_res;

  void SetUp() override {
    const size_t vect_sz = 1000000;
    input_data.resize(vect_sz);

    for (size_t i = 0; i < vect_sz; ++i) {
      auto r = (i * 1103515245ULL + 12345ULL) % 1000000;
      input_data[i] = (static_cast<double>(r) / 100.0) - 5000.0;
    }

    expected_res = input_data;
    std::ranges::sort(expected_res);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_res.size()) {
      return false;
    }
    return output_data == expected_res;
  }
};

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, FrolovaSRadixSortDoubleSEQ>(PPC_SETTINGS_example_threads);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = FrolovaSRadixSortDoubleRunPerfTests::CustomPerfTestName;

TEST_P(FrolovaSRadixSortDoubleRunPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

INSTANTIATE_TEST_SUITE_P(RunPerf, FrolovaSRadixSortDoubleRunPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace frolova_s_radix_sort_double
