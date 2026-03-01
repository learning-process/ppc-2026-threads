#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"
#include "zhurin_i_gauss_kernel_seq/common/include/common.hpp"
#include "zhurin_i_gauss_kernel_seq/seq/include/ops_seq.hpp"

namespace zhurin_i_gauss_kernel_seq {

struct FuncTestCase {
  int id;
  InType input;
  OutType expected;
};

class ZhurinIGaussKernelFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, FuncTestCase> {
 public:
  static std::string PrintTestName(
      const testing::TestParamInfo<std::tuple<std::function<std::shared_ptr<ppc::task::Task<InType, OutType>>(InType)>,
                                              std::string, FuncTestCase>> &info) {
    const auto &task_name = std::get<1>(info.param);
    const auto &test_param = std::get<2>(info.param);
    return task_name + "_Test" + std::to_string(test_param.id);
  }

 protected:
  void SetUp() override {
    const auto &params = std::get<2>(GetParam());
    input_data_ = params.input;
    expected_output_ = params.expected;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_output_.size()) {
      return false;
    }
    for (size_t i = 0; i < output_data.size(); ++i) {
      if (output_data[i].size() != expected_output_[i].size()) {
        return false;
      }
      for (size_t j = 0; j < output_data[i].size(); ++j) {
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

InType make_input(int w, int h, int p, std::vector<std::vector<int>> img) {
  return std::make_tuple(w, h, p, img);
}

const std::array<FuncTestCase, 6> kAllTests = {
    {{1, make_input(1, 1, 1, {{16}}), {{4}}},
     {2,
      make_input(3, 3, 1, std::vector<std::vector<int>>(3, std::vector<int>(3, 1))),
      {{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}},
     {3, make_input(3, 3, 1, {{0, 0, 0}, {0, 16, 0}, {0, 0, 0}}), {{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}},
     {4, make_input(3, 3, 2, {{16, 0, 0}, {0, 0, 0}, {0, 0, 0}}), {{4, 2, 0}, {2, 1, 0}, {0, 0, 0}}},
     {5, make_input(2, 2, 1, {{1, 2}, {3, 4}}), {{1, 1}, {1, 1}}},
     {6, make_input(4, 4, 4, std::vector<std::vector<int>>(4, std::vector<int>(4, 0))),
      OutType(4, std::vector<int>(4, 0))}}};

const auto kAllTasksList =
    ppc::util::AddFuncTask<ZhurinIGaussKernelSEQ, InType>(kAllTests, PPC_SETTINGS_zhurin_i_gauss_kernel_seq);

inline const auto kGtestValues = ppc::util::ExpandToValues(kAllTasksList);

TEST_P(ZhurinIGaussKernelFuncTests, AllTests) {
  ExecuteTest(GetParam());
}

INSTANTIATE_TEST_SUITE_P(ZhurinIGaussKernel, ZhurinIGaussKernelFuncTests, kGtestValues,
                         ZhurinIGaussKernelFuncTests::PrintTestName);

TEST(ZhurinIGaussKernelNegativeTest, InvalidWidth) {
  int width = 0, height = 3, parts = 1;
  std::vector<std::vector<int>> img(height, std::vector<int>(3, 0));
  InType in = std::make_tuple(width, height, parts, img);
  auto task = std::make_shared<ZhurinIGaussKernelSEQ>(in);
  EXPECT_FALSE(task->Validation());
}

TEST(ZhurinIGaussKernelNegativeTest, InvalidHeight) {
  int width = 3, height = -1, parts = 1;
  std::vector<std::vector<int>> img(1, std::vector<int>(3, 0));
  InType in = std::make_tuple(width, height, parts, img);
  auto task = std::make_shared<ZhurinIGaussKernelSEQ>(in);
  EXPECT_FALSE(task->Validation());
}

TEST(ZhurinIGaussKernelNegativeTest, InvalidPartsZero) {
  int width = 3, height = 3, parts = 0;
  std::vector<std::vector<int>> img(height, std::vector<int>(width, 0));
  InType in = std::make_tuple(width, height, parts, img);
  auto task = std::make_shared<ZhurinIGaussKernelSEQ>(in);
  EXPECT_FALSE(task->Validation());
}

TEST(ZhurinIGaussKernelNegativeTest, InvalidPartsTooLarge) {
  int width = 3, height = 3, parts = 5;
  std::vector<std::vector<int>> img(height, std::vector<int>(width, 0));
  InType in = std::make_tuple(width, height, parts, img);
  auto task = std::make_shared<ZhurinIGaussKernelSEQ>(in);
  EXPECT_FALSE(task->Validation());
}

TEST(ZhurinIGaussKernelNegativeTest, ImageRowsMismatch) {
  int width = 3, height = 3, parts = 1;
  std::vector<std::vector<int>> img(2, std::vector<int>(width, 0));
  InType in = std::make_tuple(width, height, parts, img);
  auto task = std::make_shared<ZhurinIGaussKernelSEQ>(in);
  EXPECT_FALSE(task->Validation());
}

TEST(ZhurinIGaussKernelNegativeTest, ImageColsMismatch) {
  int width = 3, height = 3, parts = 1;
  std::vector<std::vector<int>> img(height, std::vector<int>(2, 0));
  InType in = std::make_tuple(width, height, parts, img);
  auto task = std::make_shared<ZhurinIGaussKernelSEQ>(in);
  EXPECT_FALSE(task->Validation());
}

}  // namespace

}  // namespace zhurin_i_gauss_kernel_seq
