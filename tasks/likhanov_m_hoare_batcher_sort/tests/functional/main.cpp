#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <numeric>
#include <random>
#include <string>
#include <tuple>

#include "likhanov_m_hoare_batcher_sort/common/include/common.hpp"
#include "likhanov_m_hoare_batcher_sort/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace likhanov_m_hoare_batcher_sort {

class LikhanovMHoareBatcherSortRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    const auto &full_param = GetParam();
    const TestType &test_param = std::get<2>(full_param);

    const int size = std::get<0>(test_param);
    const std::string type = std::get<1>(test_param);

    input_data_.resize(size);

    std::ranges::iota(input_data_, 0);

    if (type == "reverse") {
      std::ranges::reverse(input_data_);

    } else if (type == "random") {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_int_distribution<int> dist(-1000, 1000);

      for (auto &v : input_data_) {
        v = dist(gen);
      }

    } else if (type == "duplicates") {
      for (int i = 0; i < size; ++i) {
        input_data_[i] = i % 5;
      }

    } else if (type == "negative") {
      for (int i = 0; i < size; ++i) {
        input_data_[i] = -i;
      }
    }

    expected_data_ = input_data_;
    std::ranges::sort(expected_data_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == expected_data_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_data_;
};

namespace {

const std::array<TestType, 15> kTestParam = {
    std::make_tuple(1, "single"),   std::make_tuple(2, "two"),         std::make_tuple(3, "three"),
    std::make_tuple(5, "reverse"),  std::make_tuple(10, "sorted"),     std::make_tuple(10, "reverse"),
    std::make_tuple(10, "random"),  std::make_tuple(10, "duplicates"), std::make_tuple(10, "negative"),
    std::make_tuple(100, "random"), std::make_tuple(101, "random"),    std::make_tuple(256, "reverse"),
    std::make_tuple(512, "random"), std::make_tuple(1000, "random"),   std::make_tuple(2048, "random")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<LikhanovMHoareBatcherSortSEQ, InType>(kTestParam, PPC_SETTINGS_example_threads));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    LikhanovMHoareBatcherSortRunFuncTests::PrintFuncTestName<LikhanovMHoareBatcherSortRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(HoareBatcherFuncTests, LikhanovMHoareBatcherSortRunFuncTests, kGtestValues, kPerfTestName);

TEST(HoareBatcherSortSEQ, EmptyVectorValidation) {
  InType input;
  LikhanovMHoareBatcherSortSEQ task(input);
  ASSERT_FALSE(task.Validation());
}

TEST(HoareBatcherSortSEQ, AlreadySortedManual) {
  InType input{1, 2, 3, 4, 5};
  LikhanovMHoareBatcherSortSEQ task(input);

  ASSERT_TRUE(task.Validation());
  task.PreProcessing();
  task.Run();
  task.PostProcessing();

  EXPECT_EQ(task.GetOutput(), input);
}

TEST(HoareBatcherSortSEQ, AllElementsEqual) {
  InType input(100, 42);
  LikhanovMHoareBatcherSortSEQ task(input);

  ASSERT_TRUE(task.Validation());
  task.PreProcessing();
  task.Run();
  task.PostProcessing();

  EXPECT_EQ(task.GetOutput(), input);
}

TEST(HoareBatcherSortSEQ, LargeRandomStressTest) {
  const int size = 10000;

  InType input(size);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dist(-100000, 100000);

  for (auto &v : input) {
    v = dist(gen);
  }

  InType expected = input;
  std::ranges::sort(expected);

  LikhanovMHoareBatcherSortSEQ task(input);

  ASSERT_TRUE(task.Validation());
  task.PreProcessing();
  task.Run();
  task.PostProcessing();

  EXPECT_EQ(task.GetOutput(), expected);
}

TEST_P(LikhanovMHoareBatcherSortRunFuncTests, RunFunctionalTest) {
  ExecuteTest(GetParam());
}

}  // namespace
}  // namespace likhanov_m_hoare_batcher_sort
