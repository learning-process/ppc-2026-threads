#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <random>
#include <string>
#include <tuple>

#include "nikitina_v_hoar_sort_batcher/common/include/common.hpp"
#include "nikitina_v_hoar_sort_batcher/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace nikitina_v_hoar_sort_batcher {

class NikitinaVHoarSortBatcherFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int test_id = std::get<0>(params);
    input_data_.clear();
    if (test_id == 1) {
      input_data_ = {};
    } else if (test_id == 2) {
      input_data_ = {1, 2, 3, 4, 5};
    } else if (test_id == 3) {
      input_data_ = {5, 4, 3, 2, 1};
    } else if (test_id == 4) {
      std::mt19937 gen(42);
      std::uniform_int_distribution<> dist(-100, 100);
      input_data_.resize(20);
      for (int &x : input_data_) {
        x = dist(gen);
      }
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

namespace {

TEST_P(NikitinaVHoarSortBatcherFuncTests, RunFuncTests) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 4> kTestParam = {std::make_tuple(1, "empty"), std::make_tuple(2, "sorted"),
                                            std::make_tuple(3, "reverse"), std::make_tuple(4, "random")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<HoareSortBatcherSEQ, InType>(kTestParam, PPC_SETTINGS_nikitina_v_hoar_sort_batcher));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = NikitinaVHoarSortBatcherFuncTests::PrintFuncTestName<NikitinaVHoarSortBatcherFuncTests>;

INSTANTIATE_TEST_SUITE_P(NikitinaVHoarSortBatcherTests, NikitinaVHoarSortBatcherFuncTests, kGtestValues, kTestName);

}  // namespace
}  // namespace nikitina_v_hoar_sort_batcher
