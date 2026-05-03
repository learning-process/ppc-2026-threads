#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <random>
#include <string>
#include <tuple>

#include "shekhirev_v_hoare_batcher_sort_seq/common/include/common.hpp"
#include "shekhirev_v_hoare_batcher_sort_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace shekhirev_v_hoare_batcher_sort_seq {

using TaskTestType = std::tuple<size_t, int>;

class ShekhirevVFuncTest : public ppc::util::BaseRunFuncTests<InType, OutType, TaskTestType> {
 public:
  static std::string PrintTestParam(const TaskTestType &test_param) {
    auto [size, seed] = test_param;
    return "Size_" + std::to_string(size) + "_Seed_" + std::to_string(seed);
  }

 protected:
  void SetUp() override {
    auto [size, seed] = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_.resize(size);
    if (size > 0) {
      std::mt19937 gen(seed);
      std::uniform_int_distribution<int> dist(-1000, 1000);
      for (auto &val : input_data_) {
        val = dist(gen);
      }
    }
  }

  InType GetTestInputData() final {
    return input_data_;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::ranges::is_sorted(output_data);
  }

 private:
  InType input_data_;
};

TEST_P(ShekhirevVFuncTest, SeqSortTests) {
  ExecuteTest(GetParam());
}

namespace {
const std::array<TaskTestType, 6> kTestParams = {std::make_tuple(0, 42),  std::make_tuple(1, 42),
                                                 std::make_tuple(8, 7),   std::make_tuple(13, 13),
                                                 std::make_tuple(128, 1), std::make_tuple(200, 123)};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<ShekhirevHoareBatcherSortSEQ, InType>(
    kTestParams, PPC_SETTINGS_shekhirev_v_hoare_batcher_sort_seq));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

INSTANTIATE_TEST_SUITE_P(SeqSortTests_Group, ShekhirevVFuncTest, kGtestValues,
                         ShekhirevVFuncTest::PrintFuncTestName<ShekhirevVFuncTest>);
}  // namespace

}  // namespace shekhirev_v_hoare_batcher_sort_seq
