#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <random>
#include <string>
#include <vector>

#include "zaharov_g_linear_contrast_stretch/common/include/common.hpp"
#include "zaharov_g_linear_contrast_stretch/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

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

class ZaharovGRunPerfTestsLinContrStr : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    constexpr size_t kSize = 1u << 20;
    input_data_.resize(kSize);

    std::mt19937 gen(20260301u);
    std::uniform_int_distribution<int> dist(0, 255);
    for (auto& v : input_data_) {
      v = static_cast<uint8_t>(dist(gen));
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

TEST_P(ZaharovGRunPerfTestsLinContrStr, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ZaharovGLinContrStrSEQ>(GetSettingsPath());

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = ZaharovGRunPerfTestsLinContrStr::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ZaharovGRunPerfTestsLinContrStr, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace zaharov_g_linear_contrast_stretch
