// func_tests.cpp
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <string>
#include <tuple>
#include <vector>

#include "redkina_a_sort_hoar_batcher_seq/common/include/common.hpp"
#include "redkina_a_sort_hoar_batcher_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace redkina_a_sort_hoar_batcher_seq {

class RedkinaASortHoarBatcherFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_size" + std::to_string(std::get<1>(test_param).size());
  }

 protected:
  void SetUp() override {
    auto params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<1>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != input_data_.size()) {
      return false;
    }
    if (!std::ranges::is_sorted(output_data)) {
      return false;
    }

    std::vector<int> input_copy = input_data_;
    std::vector<int> output_copy = output_data;
    std::ranges::sort(input_copy);
    std::ranges::sort(output_copy);
    return input_copy == output_copy;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

std::vector<int> MakeVector(std::initializer_list<int> list) {
  return {list};
}

const std::array<TestType, 23> kTestCases = {{std::make_tuple(1, MakeVector({})),
                                              std::make_tuple(2, MakeVector({0})),
                                              std::make_tuple(3, MakeVector({7})),
                                              std::make_tuple(4, MakeVector({0, 1})),
                                              std::make_tuple(5, MakeVector({7, 8})),
                                              std::make_tuple(6, MakeVector({1, 0})),
                                              std::make_tuple(7, MakeVector({9, 7})),
                                              std::make_tuple(8, MakeVector({0, 1, 7})),
                                              std::make_tuple(9, MakeVector({0, 7, 1})),
                                              std::make_tuple(10, MakeVector({1, 0, 7})),
                                              std::make_tuple(11, MakeVector({1, 7, 0})),
                                              std::make_tuple(12, MakeVector({7, 0, 1})),
                                              std::make_tuple(13, MakeVector({7, 1, 0})),
                                              std::make_tuple(14, MakeVector({7, 8, 9})),
                                              std::make_tuple(15, MakeVector({9, 8, 7})),
                                              std::make_tuple(16, MakeVector({7, 9, 8})),
                                              std::make_tuple(17, MakeVector({0, 0, 1, 1, 7, 7})),
                                              std::make_tuple(18, MakeVector({7, 7, 7, 7})),
                                              std::make_tuple(19, MakeVector({8, 8, 7, 7, 9, 9})),
                                              std::make_tuple(20, MakeVector({0, 1, 7, 8, 9, 0, 1, 7, 8, 9})),
                                              std::make_tuple(21, MakeVector({9, 8, 7, 1, 0, 9, 8, 7, 1, 0})),
                                              std::make_tuple(22, MakeVector({7, 8, 9, 0, 1, 7, 8, 9, 0, 1})),
                                              std::make_tuple(23, [] {
  std::vector<int> v;
  v.reserve(20);
  for (int i = 0; i < 20; ++i) {
    int remainder = i % 5;
    if (remainder == 0) {
      v.push_back(0);
    } else if (remainder == 1) {
      v.push_back(1);
    } else if (remainder == 2) {
      v.push_back(7);
    } else if (remainder == 3) {
      v.push_back(8);
    } else {
      v.push_back(9);
    }
  }
  return v;
}())}};

const auto kTestTasksList = ppc::util::AddFuncTask<RedkinaASortHoarBatcherSEQ, InType>(
    kTestCases, PPC_SETTINGS_redkina_a_sort_hoar_batcher_seq);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = RedkinaASortHoarBatcherFuncTests::PrintFuncTestName<RedkinaASortHoarBatcherFuncTests>;

INSTANTIATE_TEST_SUITE_P(SortHoarBatcherTests, RedkinaASortHoarBatcherFuncTests, kGtestValues, kTestName);

TEST_P(RedkinaASortHoarBatcherFuncTests, SortCheck) {
  ExecuteTest(GetParam());
}

}  // namespace

}  // namespace redkina_a_sort_hoar_batcher_seq
