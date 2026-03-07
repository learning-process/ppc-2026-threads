#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "zaharov_g_linear_contrast_stretch/common/include/common.hpp"
#include "zaharov_g_linear_contrast_stretch/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace zaharov_g_linear_contrast_stretch {

namespace {

std::string GetSettingsPath() {
  return (std::filesystem::path(PPC_PATH_TO_PROJECT) / "tasks" / "zaharov_g_linear_contrast_stretch" /
          "settings.json")
      .string();
}

OutType ReferenceLinContrStr(const InType& input) {
  OutType output;
  if (input.empty()) {
    return output;
  }

  output.resize(input.size());
  auto [min_it, max_it] = std::minmax_element(input.begin(), input.end());
  const uint8_t min_el = *min_it;
  const uint8_t max_el = *max_it;

  if (max_el > min_el) {
    const int denom = static_cast<int>(max_el) - static_cast<int>(min_el);
    for (size_t i = 0; i < input.size(); ++i) {
      const int value = (static_cast<int>(input[i]) - static_cast<int>(min_el)) * 255 / denom;
      output[i] = static_cast<uint8_t>(std::clamp(value, 0, 255));
    }
  } else {
    output.assign(input.begin(), input.end());
  }
  return output;
}

}  // namespace

TEST(ZaharovGLinContrStrSEQ, ValidationRejectsEmptyInput) {
  const InType input;
  ZaharovGLinContrStrSEQ task(input);
  EXPECT_FALSE(task.Validation());
}

class ZaharovGRunFuncTestsLinContrStr : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType& test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    auto test_param = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    const int n = std::get<0>(test_param);
    const std::string& pattern = std::get<1>(test_param);

    if (pattern == "shifted") {
      input_data_ = {50, 75, 100};
    } else if (pattern == "constant") {
      input_data_.assign(static_cast<size_t>(n), 7);
    } else if (pattern == "full_range") {
      input_data_.resize(256);
      for (size_t i = 0; i < input_data_.size(); ++i) {
        input_data_[i] = static_cast<uint8_t>(i);
      }
    } else if (pattern == "random") {
      input_data_.resize(static_cast<size_t>(n));
      std::mt19937 gen(12345u + static_cast<uint32_t>(n));
      std::uniform_int_distribution<int> dist(0, 255);
      for (auto& v : input_data_) {
        v = static_cast<uint8_t>(dist(gen));
      }
    } else {
      FAIL() << "Unknown pattern: " << pattern;
    }

    expected_data_ = ReferenceLinContrStr(input_data_);
  }

  bool CheckTestOutputData(OutType& output_data) final {
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

TEST_P(ZaharovGRunFuncTestsLinContrStr, Correctness) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 4> kCases = {
    std::make_tuple(3, "shifted"),
    std::make_tuple(10, "constant"),
    std::make_tuple(256, "full_range"),
    std::make_tuple(1000, "random"),
};

const auto kTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<ZaharovGLinContrStrSEQ, InType>(kCases, GetSettingsPath()));

const auto kGtestValues = ppc::util::ExpandToValues(kTasksList);

const auto kTestName = ZaharovGRunFuncTestsLinContrStr::PrintFuncTestName<ZaharovGRunFuncTestsLinContrStr>;

INSTANTIATE_TEST_SUITE_P(LinContrStrTests, ZaharovGRunFuncTestsLinContrStr, kGtestValues, kTestName);

}  // namespace

}  // namespace zaharov_g_linear_contrast_stretch
