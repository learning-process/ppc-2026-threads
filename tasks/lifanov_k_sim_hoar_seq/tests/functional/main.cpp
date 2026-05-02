#include <gtest/gtest.h>

#include <algorithm>
#include <random>
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

static std::vector<int> GetRandomVector(size_t size) {
  std::vector<int> vec(size);
  std::random_device rd;
  std::mt19937 gen(rd());
  for (size_t i = 0; i < size; ++i) {
    vec[i] = static_cast<int>(gen() % 1000);
  }
  return vec;
}

static std::vector<int> GetSortedVector(std::vector<int> vec) {
  std::sort(vec.begin(), vec.end());
  return vec;
}

const std::vector<int> v1 = GetRandomVector(10);
const std::vector<int> v2 = GetRandomVector(100);
const std::vector<int> v3 = {5, 4, 3, 2, 1};

const std::array<TestType, 3> kTestParam = {std::make_tuple(v1, GetSortedVector(v1), "size_10"),
                                            std::make_tuple(v2, GetSortedVector(v2), "size_100"),
                                            std::make_tuple(v3, GetSortedVector(v3), "reverse_5")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<LifanovKSimpleHoarSEQ, InType>(kTestParam, "lifanov_k_sim_hoar_seq");

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = LifanovKRunFuncTestsSEQ::PrintFuncTestName<LifanovKRunFuncTestsSEQ>;

INSTANTIATE_TEST_SUITE_P(HoarSortTests, LifanovKRunFuncTestsSEQ, kGtestValues, kTestName);

}  // namespace
}  // namespace lifanov_k_sim_hoar_seq
