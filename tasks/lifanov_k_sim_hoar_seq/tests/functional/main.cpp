#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "lifanov_k_sim_hoar_seq/common/include/common.hpp"
#include "lifanov_k_sim_hoar_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace lifanov_k_sim_hoar_seq {

class LifanovKRunFuncTestsSEQ : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    auto params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_output_ = std::get<1>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return expected_output_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_output_;
};

namespace {

TEST_P(LifanovKRunFuncTestsSEQ, SortTest) {
  ExecuteTest(GetParam());
}

std::vector<int> GetRandomVector(std::size_t size) {
  std::vector<int> vec(size);
  std::random_device rd;
  std::mt19937 gen(rd());

  for (std::size_t i = 0; i < size; ++i) {
    vec[i] = static_cast<int>(gen() % 1000);
  }

  return vec;
}

std::vector<int> GetSortedVector(std::vector<int> vec) {
  std::ranges::sort(vec);
  return vec;
}

const std::vector<int> kV1 = GetRandomVector(10);
const std::vector<int> kV2 = GetRandomVector(100);
const std::vector<int> kV3 = {5, 4, 3, 2, 1};

const std::array<TestType, 3> kTestParam = {std::make_tuple(kV1, GetSortedVector(kV1), "size_10"),
                                            std::make_tuple(kV2, GetSortedVector(kV2), "size_100"),
                                            std::make_tuple(kV3, GetSortedVector(kV3), "reverse_5")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<LifanovKSimpleHoarSEQ, InType>(kTestParam, PPC_SETTINGS_lifanov_k_sim_hoar_seq);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = LifanovKRunFuncTestsSEQ::PrintFuncTestName<LifanovKRunFuncTestsSEQ>;

INSTANTIATE_TEST_SUITE_P(HoarSortTests, LifanovKRunFuncTestsSEQ, kGtestValues, kTestName);

}  // namespace
}  // namespace lifanov_k_sim_hoar_seq
