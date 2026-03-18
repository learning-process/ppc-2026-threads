#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <random>
#include <tuple>
#include <vector>

#include "romanov_a_gauss_block/common/include/common.hpp"
#include "romanov_a_gauss_block/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace romanov_a_gauss_block {

class RomanovAPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kWidth_ = 7680;
  const int kHeight_ = 4320;
  InType input_data_;

  void SetUp() override {
    std::vector<uint8_t> picture(static_cast<size_t>(kWidth_ * kHeight_ * 3));
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> dist(0, 255);
    for (uint8_t &v : picture) {
      v = static_cast<uint8_t>(dist(rng));
    }
    input_data_ = std::make_tuple(kWidth_, kHeight_, picture);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::get<2>(input_data_).size() == output_data.size();
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(RomanovAPerfTestThreads, GaussFilter) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, RomanovAGaussBlockSEQ>(PPC_SETTINGS_example_threads);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = RomanovAPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(GaussFilter, RomanovAPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace romanov_a_gauss_block
