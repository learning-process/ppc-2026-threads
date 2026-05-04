#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <iostream>
#include <string>
#include <tuple>

#include "artyushkina_markirovka/common/include/common.hpp"
#include "artyushkina_markirovka/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace artyushkina_markirovka {

class ArtyushkinaMarkirovkaFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int test_id = std::get<0>(params);

    switch (test_id) {
      case 0: {
        input_data_ = {3, 3, 0, 0, 255, 0, 255, 255, 255, 255, 0};
        expected_ = {3, 3, 1, 1, 0, 1, 0, 0, 0, 0, 2};
        break;
      }
      case 1: {
        input_data_ = {3, 3, 0, 0, 0, 0, 255, 0, 0, 0, 255};
        expected_ = {3, 3, 1, 1, 1, 1, 0, 1, 1, 1, 0};
        break;
      }
      case 2: {
        input_data_ = {2, 3, 255, 255, 255, 255, 255, 255};
        expected_ = {2, 3, 0, 0, 0, 0, 0, 0};
        break;
      }
      case 3: {
        input_data_ = {2, 2, 0, 0, 0, 0};
        expected_ = {2, 2, 1, 1, 1, 1};
        break;
      }
      case 4: {
        input_data_ = {3, 4, 0, 0, 0, 0, 255, 255, 255, 255, 0, 0, 0, 0};
        expected_ = {3, 4, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2};
        break;
      }
      default:
        break;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data != expected_) {
      std::cout << "Expected: ";
      for (auto val : expected_) {
        std::cout << static_cast<int>(val) << ' ';
      }
      std::cout << '\n';

      std::cout << "Actual  : ";
      for (auto val : output_data) {
        std::cout << val << ' ';
      }
      std::cout << '\n';
    }

    return output_data == expected_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_;
};

namespace {

TEST_P(ArtyushkinaMarkirovkaFuncTests, MarkingComponents) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParam = {
    std::make_tuple(0, "L_shaped_component_8connectivity"), std::make_tuple(1, "diagonal_connected_components"),
    std::make_tuple(2, "all_background"), std::make_tuple(3, "all_objects"), std::make_tuple(4, "two_horizontal_bars")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<MarkingComponentsSEQ, InType>(kTestParam, PPC_SETTINGS_artyushkina_markirovka);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = ArtyushkinaMarkirovkaFuncTests::PrintFuncTestName<ArtyushkinaMarkirovkaFuncTests>;

INSTANTIATE_TEST_SUITE_P(ComponentLabeling, ArtyushkinaMarkirovkaFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace artyushkina_markirovka
