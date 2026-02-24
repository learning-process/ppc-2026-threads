#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "kamaletdinov_r_bitwise_int_seq/common/include/common.hpp"
#include "kamaletdinov_r_bitwise_int_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kamaletdinov_r_bitwise_int_seq {

class KamaletdinovRBitwiseIntFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (input_data_ == output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = 0;
};

namespace {

void RunPipelineSteps(const std::shared_ptr<KamaletdinovRBitwiseIntSEQ> &task) {
  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
}

void RunTaskPipeline(const std::shared_ptr<KamaletdinovRBitwiseIntSEQ> &task, int expected) {
  RunPipelineSteps(task);
  EXPECT_TRUE(task->PostProcessing());
  EXPECT_EQ(task->GetOutput(), expected);
}

TEST_P(KamaletdinovRBitwiseIntFuncTests, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(3, "3"), std::make_tuple(5, "5"), std::make_tuple(7, "7")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<KamaletdinovRBitwiseIntSEQ, InType>(
    kTestParam, PPC_SETTINGS_kamaletdinov_r_bitwise_int_seq));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = KamaletdinovRBitwiseIntFuncTests::PrintFuncTestName<KamaletdinovRBitwiseIntFuncTests>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, KamaletdinovRBitwiseIntFuncTests, kGtestValues, kPerfTestName);

TEST(KamaletdinovRBitwiseIntEdge, PipelineEmptyInput) {
  RunTaskPipeline(std::make_shared<KamaletdinovRBitwiseIntSEQ>(0), 0);
}

TEST(KamaletdinovRBitwiseIntEdge, PipelineSingleElement) {
  RunTaskPipeline(std::make_shared<KamaletdinovRBitwiseIntSEQ>(1), 1);
}

TEST(KamaletdinovRBitwiseIntEdge, PipelineTwoElements) {
  RunTaskPipeline(std::make_shared<KamaletdinovRBitwiseIntSEQ>(2), 2);
}

TEST(KamaletdinovRBitwiseIntEdge, PipelineLargerInput) {
  RunTaskPipeline(std::make_shared<KamaletdinovRBitwiseIntSEQ>(100), 100);
}

TEST(KamaletdinovRBitwiseIntEdge, SortEmptyVector) {
  std::vector<int> data;
  BitwiseSort(data);
  EXPECT_TRUE(data.empty());
}

TEST(KamaletdinovRBitwiseIntEdge, SortSingleElement) {
  std::vector<int> data = {42};
  BitwiseSort(data);
  ASSERT_EQ(data.size(), 1U);
  EXPECT_EQ(data[0], 42);
}

TEST(KamaletdinovRBitwiseIntEdge, SortAlreadySorted) {
  std::vector<int> data = {1, 2, 3, 4, 5};
  BitwiseSort(data);
  EXPECT_TRUE(std::ranges::is_sorted(data));
  std::vector<int> expected = {1, 2, 3, 4, 5};
  EXPECT_EQ(data, expected);
}

TEST(KamaletdinovRBitwiseIntEdge, SortReverseSorted) {
  std::vector<int> data = {5, 4, 3, 2, 1};
  BitwiseSort(data);
  std::vector<int> expected = {1, 2, 3, 4, 5};
  EXPECT_EQ(data, expected);
}

TEST(KamaletdinovRBitwiseIntEdge, SortAllSameElements) {
  std::vector<int> data = {7, 7, 7, 7};
  BitwiseSort(data);
  std::vector<int> expected = {7, 7, 7, 7};
  EXPECT_EQ(data, expected);
}

TEST(KamaletdinovRBitwiseIntEdge, SortAllNegative) {
  std::vector<int> data = {-3, -1, -5, -2, -4};
  BitwiseSort(data);
  std::vector<int> expected = {-5, -4, -3, -2, -1};
  EXPECT_EQ(data, expected);
}

TEST(KamaletdinovRBitwiseIntEdge, SortMixedPositiveNegative) {
  std::vector<int> data = {3, -1, 4, -1, 5, -9, 2, 6};
  BitwiseSort(data);
  std::vector<int> expected = {-9, -1, -1, 2, 3, 4, 5, 6};
  EXPECT_EQ(data, expected);
}

TEST(KamaletdinovRBitwiseIntEdge, SortWithZeros) {
  std::vector<int> data = {0, 3, 0, -2, 0, 1};
  BitwiseSort(data);
  std::vector<int> expected = {-2, 0, 0, 0, 1, 3};
  EXPECT_EQ(data, expected);
}

TEST(KamaletdinovRBitwiseIntEdge, SortTwoElements) {
  std::vector<int> data = {2, 1};
  BitwiseSort(data);
  std::vector<int> expected = {1, 2};
  EXPECT_EQ(data, expected);
}

TEST(KamaletdinovRBitwiseIntEdge, SortLargeNumbers) {
  std::vector<int> data = {1000000, 999999, 1, 500000};
  BitwiseSort(data);
  std::vector<int> expected = {1, 500000, 999999, 1000000};
  EXPECT_EQ(data, expected);
}

TEST(KamaletdinovRBitwiseIntEdge, SortDuplicatesWithNegatives) {
  std::vector<int> data = {-3, 5, -3, 5, 0, 0};
  BitwiseSort(data);
  std::vector<int> expected = {-3, -3, 0, 0, 5, 5};
  EXPECT_EQ(data, expected);
}

TEST(KamaletdinovRBitwiseIntEdge, SortAllZeros) {
  std::vector<int> data = {0, 0, 0, 0};
  BitwiseSort(data);
  std::vector<int> expected = {0, 0, 0, 0};
  EXPECT_EQ(data, expected);
}

TEST(KamaletdinovRBitwiseIntEdge, SortTwoNegative) {
  std::vector<int> data = {-1, -2};
  BitwiseSort(data);
  std::vector<int> expected = {-2, -1};
  EXPECT_EQ(data, expected);
}

}  // namespace

}  // namespace kamaletdinov_r_bitwise_int_seq
