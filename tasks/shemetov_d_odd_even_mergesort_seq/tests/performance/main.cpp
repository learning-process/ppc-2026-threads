#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <random>
#include <vector>

#include "shemetov_d_odd_even_mergesort_seq/common/include/common.hpp"
#include "shemetov_d_odd_even_mergesort_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace shemetov_d_odd_even_mergesort {

class ShemetovDRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    const size_t size = (1ULL << 17) + 1;
    std::vector<int> gen_array(size);

    std::random_device rnd_dvc;
    std::mt19937 gen(rnd_dvc());
    std::uniform_int_distribution<int> dist(-1000000, 1000000);

    std::ranges::generate(gen_array.begin(), gen_array.end(), [&]() { return dist(gen); });

    test_array_ = {size, gen_array};
  }

  InType GetTestInputData() final {
    return test_array_;
  }

  bool CheckTestOutputData(OutType &output) final {
    return std::ranges::is_sorted(output.begin(), output.end());
  }

 private:
  InType test_array_;
};

TEST_P(ShemetovDRunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ShemetovDOddEvenMergeSortSEQ>(PPC_SETTINGS_shemetov_d_odd_even_mergesort_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ShemetovDRunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(OddEvenMergeSortPerfTest, ShemetovDRunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace shemetov_d_odd_even_mergesort
