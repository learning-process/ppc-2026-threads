#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <iostream>
#include <string>
#include <tuple>

#include "artyushkina_markirovka/common/include/common.hpp"
#include "artyushkina_markirovka/omp/include/ops_omp.hpp"
#include "artyushkina_markirovka/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace artyushkina_markirovka {
namespace {

void PrintMatrix(const OutType &data, const std::string &title) {
  std::cout << title;
  int rows = static_cast<int>(data[0]);
  int cols = static_cast<int>(data[1]);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      std::size_t idx =
          (static_cast<std::size_t>(i) * static_cast<std::size_t>(cols)) + static_cast<std::size_t>(j) + 2;
      std::cout << data[idx] << " ";
    }
    std::cout << "\n";
  }
  std::cout << "=====================\n\n";
}

void PrintInputMatrix(const InType &data) {
  int rows = static_cast<int>(data[0]);
  int cols = static_cast<int>(data[1]);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      std::size_t idx =
          (static_cast<std::size_t>(i) * static_cast<std::size_t>(cols)) + static_cast<std::size_t>(j) + 2;
      std::cout << static_cast<int>(data[idx]) << " ";
    }
    std::cout << "\n";
  }
}

}  // namespace

class ArtyushkinaMarkirovkaFuncTestsSEQ : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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

class ArtyushkinaMarkirovkaFuncTestsOMP : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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
      case 5: {
        input_data_ = {4, 4, 0, 0, 255, 255, 0, 255, 0, 255, 255, 0, 0, 255, 0, 0, 255, 0};
        expected_ = {4, 4, 1, 1, 0, 0, 1, 0, 2, 0, 0, 2, 2, 0, 3, 3, 0, 4};
        break;
      }
      case 6: {
        input_data_ = {2, 2, 0, 255, 255, 0};
        expected_ = {2, 2, 1, 0, 0, 1};
        break;
      }
      default:
        break;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    int test_id = std::get<0>(std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam()));

    if (test_id == 5 || test_id == 6) {
      std::cout << "\n=== INPUT MATRIX (test " << test_id << ") ===\n";
      PrintInputMatrix(input_data_);
      std::cout << "========================\n\n";
      PrintMatrix(output_data, "=== OUTPUT LABELS ===\n");
    }

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

TEST_P(ArtyushkinaMarkirovkaFuncTestsSEQ, MarkingComponentsSEQ) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParamSEQ = {
    std::make_tuple(0, "L_shaped_component_8connectivity"), std::make_tuple(1, "diagonal_connected_components"),
    std::make_tuple(2, "all_background"), std::make_tuple(3, "all_objects"), std::make_tuple(4, "two_horizontal_bars")};

const auto kTestTasksListSEQ =
    ppc::util::AddFuncTask<MarkingComponentsSEQ, InType>(kTestParamSEQ, PPC_SETTINGS_artyushkina_markirovka);

const auto kGtestValuesSEQ = ppc::util::ExpandToValues(kTestTasksListSEQ);

const auto kPerfTestNameSEQ = ArtyushkinaMarkirovkaFuncTestsSEQ::PrintFuncTestName<ArtyushkinaMarkirovkaFuncTestsSEQ>;

INSTANTIATE_TEST_SUITE_P(ComponentLabelingSEQ, ArtyushkinaMarkirovkaFuncTestsSEQ, kGtestValuesSEQ, kPerfTestNameSEQ);

TEST_P(ArtyushkinaMarkirovkaFuncTestsOMP, MarkingComponentsOMP) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 7> kTestParamOMP = {std::make_tuple(0, "L_shaped_component_8connectivity"),
                                               std::make_tuple(1, "diagonal_connected_components"),
                                               std::make_tuple(2, "all_background"),
                                               std::make_tuple(3, "all_objects"),
                                               std::make_tuple(4, "two_horizontal_bars"),
                                               std::make_tuple(5, "complex_shape_multiple_components"),
                                               std::make_tuple(6, "diagonal_connectivity_check")};

const auto kTestTasksListOMP =
    ppc::util::AddFuncTask<MarkingComponentsOMP, InType>(kTestParamOMP, PPC_SETTINGS_artyushkina_markirovka);

const auto kGtestValuesOMP = ppc::util::ExpandToValues(kTestTasksListOMP);

const auto kPerfTestNameOMP = ArtyushkinaMarkirovkaFuncTestsOMP::PrintFuncTestName<ArtyushkinaMarkirovkaFuncTestsOMP>;

INSTANTIATE_TEST_SUITE_P(ComponentLabelingOMP, ArtyushkinaMarkirovkaFuncTestsOMP, kGtestValuesOMP, kPerfTestNameOMP);

}  // namespace

}  // namespace artyushkina_markirovka
