#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "krasnopevtseva_v_hoare_batcher_sort/common/include/common.hpp"
#include "krasnopevtseva_v_hoare_batcher_sort/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace krasnopevtseva_v_hoare_batcher_sort {

class KrasnopevtsevaVRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    const auto &input = std::get<0>(test_param);
    size_t size = input.size();
    std::string s;
    for (size_t i = 0; i < size; i++) {
      s += std::to_string(input[i]);
      if (i < size - 1) {
        s += "_";
      }
    }
    std::string result = "array_" + s + "_size" + std::to_string(size) + "_" + std::get<1>(test_param);
    return result;
  }

 protected:
  void SetUp() override {
    auto test_param = std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(test_param);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    size_t size = output_data.size();
    bool result = true;
    for (size_t i = 0; i < size - 1; i++) {
      if (output_data[i] > output_data[i + 1]) {
        result = false;
      }
    }
    return result;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(KrasnopevtsevaVRunFuncTestsThreads, HoareBatcherSort) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 6> kTestParam = {
    std::make_tuple(std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}, "sorted_array"),
    std::make_tuple(std::vector<int>{5, 1, 3, 4, 2, 34, 24, 16, 31, 666, 22, 14, 52, 67, 13, 99, 9, 6, 28, 35},
                    "default_array"),
    std::make_tuple(std::vector<int>{10, 1, 30, 2}, "short_array"),
    std::make_tuple(std::vector<int>{1000, 3246, 10,   31, 120, 4,   2,    1000, 23,    34,    30,   42,    1,
                                     45,   24,   15,   32, 111, 35,  25,   252,  222,   66234, 2325, 23423, 2355,
                                     745,  579,  875,  33, 66,  345, 4666, 2490, 100,   10,    3415, 234,   22,
                                     526,  372,  8432, 21, 58,  225, 865,  23,   13333, 35,    2523, 33},
                    "long_array"),
    std::make_tuple(std::vector<int>{10, 20, 20, 10, 20, 10, 20, 10, 20}, "two_number_array"),
    std::make_tuple(std::vector<int>{10000, 9999, 9998, 4356, 3662, 3000, 2500, 2300, 2200, 2000, 1999, 1900,
                                     1843,  1000, 10,   9,    8,    7,    6,    5,    4,    3,    2,    1},
                    "end_to_begin_array")};

const auto kTestTasksList = ppc::util::AddFuncTask<KrasnopevtsevaVHoareBatcherSortSEQ, InType>(
    kTestParam, PPC_SETTINGS_krasnopevtseva_v_hoare_batcher_sort);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = KrasnopevtsevaVRunFuncTestsThreads::PrintFuncTestName<KrasnopevtsevaVRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(HoareBatcherSort, KrasnopevtsevaVRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace krasnopevtseva_v_hoare_batcher_sort
