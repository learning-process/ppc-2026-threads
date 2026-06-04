#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>

#include "egashin_k_radix_simple_merge/all/include/ops_all.hpp"
#include "egashin_k_radix_simple_merge/common/include/common.hpp"
#include "egashin_k_radix_simple_merge/omp/include/ops_omp.hpp"
#include "egashin_k_radix_simple_merge/seq/include/ops_seq.hpp"
#include "egashin_k_radix_simple_merge/stl/include/ops_stl.hpp"
#include "egashin_k_radix_simple_merge/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace egashin_k_radix_simple_merge {

class EgashinKRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    constexpr size_t kDataSize = 2'000'000;
    input_data_.resize(kDataSize);

    uint64_t state = 0x123456789ABCDEF0ULL;
    for (size_t i = 0; i < kDataSize; ++i) {
      state = (state * 2862933555777941757ULL) + 3037000493ULL;
      const auto centered = static_cast<int64_t>(state % 2'000'001ULL) - 1'000'000LL;
      input_data_[i] = static_cast<double>(centered) / 128.0;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::ranges::is_sorted(output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(EgashinKRunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, EgashinKRadixSimpleMergeSEQ, EgashinKRadixSimpleMergeOMP,
                                EgashinKRadixSimpleMergeTBB, EgashinKRadixSimpleMergeSTL, EgashinKRadixSimpleMergeALL>(
        PPC_SETTINGS_egashin_k_radix_simple_merge);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = EgashinKRunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, EgashinKRunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace egashin_k_radix_simple_merge
