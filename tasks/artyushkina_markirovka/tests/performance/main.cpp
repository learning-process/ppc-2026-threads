#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <utility>

#include "artyushkina_markirovka/common/include/common.hpp"
#include "artyushkina_markirovka/omp/include/ops_omp.hpp"
#include "artyushkina_markirovka/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace artyushkina_markirovka {

class ArtyushkinaMarkirovkaPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    const int k_size = 1000;
    input_data_.resize((static_cast<std::size_t>(k_size) * static_cast<std::size_t>(k_size)) + 2);

    input_data_[0] = static_cast<uint8_t>(k_size);
    input_data_[1] = static_cast<uint8_t>(k_size);

    for (int i = 0; i < k_size; ++i) {
      for (int j = 0; j < k_size; ++j) {
        std::size_t idx =
            (static_cast<std::size_t>(i) * static_cast<std::size_t>(k_size)) + static_cast<std::size_t>(j) + 2;
        input_data_[idx] = static_cast<uint8_t>(((i + j) % 2 == 0) ? 0 : 255);
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    int rows = static_cast<int>(output_data[0]);
    int cols = static_cast<int>(output_data[1]);

    if (std::cmp_not_equal(rows, input_data_[0]) || std::cmp_not_equal(cols, input_data_[1])) {
      std::cerr << "Size mismatch: expected " << static_cast<int>(input_data_[0]) << "x"
                << static_cast<int>(input_data_[1]) << ", got " << rows << "x" << cols << std::endl;
      return false;
    }

    bool valid = true;
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        std::size_t output_idx =
            (static_cast<std::size_t>(i) * static_cast<std::size_t>(cols)) + static_cast<std::size_t>(j) + 2;
        std::size_t input_idx =
            (static_cast<std::size_t>(i) * static_cast<std::size_t>(cols)) + static_cast<std::size_t>(j) + 2;

        if (input_data_[input_idx] == 0 && output_data[output_idx] == 0) {
          std::cerr << "Object pixel at (" << i << "," << j << ") not labeled!" << std::endl;
          valid = false;
        }
        if (input_data_[input_idx] != 0 && output_data[output_idx] != 0) {
          std::cerr << "Background pixel at (" << i << "," << j << ") incorrectly labeled!" << std::endl;
          valid = false;
        }
      }
    }
    return valid;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(ArtyushkinaMarkirovkaPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasksSEQ =
    ppc::util::MakeAllPerfTasks<InType, MarkingComponentsSEQ>(PPC_SETTINGS_artyushkina_markirovka);

const auto kGtestValuesSEQ = ppc::util::TupleToGTestValues(kAllPerfTasksSEQ);

const auto kPerfTestNameSEQ = ArtyushkinaMarkirovkaPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTestsSEQ, ArtyushkinaMarkirovkaPerfTests, kGtestValuesSEQ, kPerfTestNameSEQ);

const auto kAllPerfTasksOMP =
    ppc::util::MakeAllPerfTasks<InType, MarkingComponentsOMP>(PPC_SETTINGS_artyushkina_markirovka);

const auto kGtestValuesOMP = ppc::util::TupleToGTestValues(kAllPerfTasksOMP);

const auto kPerfTestNameOMP = ArtyushkinaMarkirovkaPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTestsOMP, ArtyushkinaMarkirovkaPerfTests, kGtestValuesOMP, kPerfTestNameOMP);

}  // namespace

}  // namespace artyushkina_markirovka
