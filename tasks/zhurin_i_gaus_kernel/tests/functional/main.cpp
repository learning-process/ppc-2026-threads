#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"
#include "util/include/func_test_util.hpp"
#include "zhurin_i_gaus_kernel/common/include/common.hpp"
#include "zhurin_i_gaus_kernel/seq/include/ops_seq.hpp"

namespace zhurin_i_gaus_kernel {

using TestCase = std::tuple<int, InType, OutType>;

class ZhurinIGausKernelFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestCase> {
 public:
  static std::string PrintTestName(
      const testing::TestParamInfo<std::tuple<std::function<std::shared_ptr<ppc::task::Task<InType, OutType>>(InType)>,
                                              std::string, TestCase>> &info) {
    const auto &task_name = std::get<1>(info.param);
    const auto &test_param = std::get<2>(info.param);
    int id = std::get<0>(test_param);
    return task_name + "_Test" + std::to_string(id);
  }

 protected:
  void SetUp() override {
    const auto &params = std::get<2>(GetParam());
    input_data_ = std::get<1>(params);
    expected_output_ = std::get<2>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_output_.size()) {
      return false;
    }
    for (std::size_t i = 0; i < output_data.size(); ++i) {
      if (output_data[i].size() != expected_output_[i].size()) {
        return false;
      }
      for (std::size_t j = 0; j < output_data[i].size(); ++j) {
        if (output_data[i][j] != expected_output_[i][j]) {
          return false;
        }
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_output_;
};

namespace {

InType MakeInput(int w, int h, int p, const std::vector<std::vector<int>> &img) {
  return std::make_tuple(w, h, p, img);
}

const std::array<TestCase, 6> kTestCases = {
    {{1, MakeInput(1, 1, 1, {{16}}), {{4}}},
     {2,
      MakeInput(3, 3, 1, std::vector<std::vector<int>>(3, std::vector<int>(3, 1))),
      {{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}},
     {3, MakeInput(3, 3, 1, {{0, 0, 0}, {0, 16, 0}, {0, 0, 0}}), {{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}},
     {4, MakeInput(3, 3, 2, {{16, 0, 0}, {0, 0, 0}, {0, 0, 0}}), {{4, 2, 0}, {2, 1, 0}, {0, 0, 0}}},
     {5, MakeInput(2, 2, 1, {{1, 2}, {3, 4}}), {{1, 1}, {1, 1}}},
     {6, MakeInput(4, 4, 4, std::vector<std::vector<int>>(4, std::vector<int>(4, 0))),
      OutType(4, std::vector<int>(4, 0))}}};

const auto kAllTasksList =
    ppc::util::AddFuncTask<ZhurinIGausKernelSEQ, InType>(kTestCases, PPC_SETTINGS_zhurin_i_gaus_kernel);

inline const auto kGtestValues = ppc::util::ExpandToValues(kAllTasksList);

TEST_P(ZhurinIGausKernelFuncTests, AllTests) {
  ExecuteTest(GetParam());
}

INSTANTIATE_TEST_SUITE_P(ZhurinIGausKernel, ZhurinIGausKernelFuncTests, kGtestValues,
                         ZhurinIGausKernelFuncTests::PrintTestName);

TEST(ZhurinIGausKernelNegativeTest, InvalidWidth) {
  int width = 0;
  int height = 3;
  int parts = 1;
  std::vector<std::vector<int>> img(height, std::vector<int>(3, 0));
  InType in = std::make_tuple(width, height, parts, img);
  auto task_seq = std::make_shared<ZhurinIGausKernelSEQ>(in);
  EXPECT_FALSE(task_seq->Validation());
}

TEST(ZhurinIGausKernelNegativeTest, InvalidHeight) {
  int width = 3;
  int height = -1;
  int parts = 1;
  std::vector<std::vector<int>> img(1, std::vector<int>(3, 0));
  InType in = std::make_tuple(width, height, parts, img);
  auto task_seq = std::make_shared<ZhurinIGausKernelSEQ>(in);
  EXPECT_FALSE(task_seq->Validation());
}

TEST(ZhurinIGausKernelNegativeTest, InvalidPartsZero) {
  int width = 3;
  int height = 3;
  int parts = 0;
  std::vector<std::vector<int>> img(height, std::vector<int>(width, 0));
  InType in = std::make_tuple(width, height, parts, img);
  auto task_seq = std::make_shared<ZhurinIGausKernelSEQ>(in);
  EXPECT_FALSE(task_seq->Validation());
}

TEST(ZhurinIGausKernelNegativeTest, InvalidPartsTooLarge) {
  int width = 3;
  int height = 3;
  int parts = 5;
  std::vector<std::vector<int>> img(height, std::vector<int>(width, 0));
  InType in = std::make_tuple(width, height, parts, img);
  auto task_seq = std::make_shared<ZhurinIGausKernelSEQ>(in);
  EXPECT_FALSE(task_seq->Validation());
}

TEST(ZhurinIGausKernelNegativeTest, ImageRowsMismatch) {
  int width = 3;
  int height = 3;
  int parts = 1;
  std::vector<std::vector<int>> img(2, std::vector<int>(width, 0));
  InType in = std::make_tuple(width, height, parts, img);
  auto task_seq = std::make_shared<ZhurinIGausKernelSEQ>(in);
  EXPECT_FALSE(task_seq->Validation());
}

TEST(ZhurinIGausKernelNegativeTest, ImageColsMismatch) {
  int width = 3;
  int height = 3;
  int parts = 1;
  std::vector<std::vector<int>> img(height, std::vector<int>(2, 0));
  InType in = std::make_tuple(width, height, parts, img);
  auto task_seq = std::make_shared<ZhurinIGausKernelSEQ>(in);
  EXPECT_FALSE(task_seq->Validation());
}

}  // namespace

}  // namespace zhurin_i_gaus_kernel
