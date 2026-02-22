#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "kapanova_s_increase_contrast_seq/common/include/common.hpp"
#include "kapanova_s_increase_contrast_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace kapanova_s_increase_contrast_seq {

class KapanovaSIncreaseContrastFuncTest : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int size = std::get<0>(params);

    input_data_.resize(static_cast<size_t>(size) * static_cast<size_t>(size));
    for (size_t i = 0; i < input_data_.size(); ++i) {
      input_data_[i] = static_cast<uint8_t>(i % 256);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (input_data_.empty() || output_data.empty()) {
      return false;
    }
    if (input_data_.size() != output_data.size()) {
      return false;
    }

    auto min_in = *std::ranges::min_element(input_data_);
    auto max_in = *std::ranges::max_element(input_data_);
    auto min_out = *std::ranges::min_element(output_data);
    auto max_out = *std::ranges::max_element(output_data);

    if (min_in == max_in) {
      return std::ranges::all_of(output_data, [min_in](uint8_t v) { return v == min_in; });
    }

    return (min_out == 0 && max_out == 255);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

class KapanovaSIncreaseContrastValidationTest : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

  void ExecuteTest(ppc::util::FuncTestParam<InType, OutType, TestType> test_param) {
    const std::string &test_name =
        std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kNameTest)>(test_param);

    this->ValidateTestName(test_name);

    const auto test_env_scope = ppc::util::test::MakePerTestEnvForCurrentGTest(test_name);

    if (this->IsTestDisabled(test_name)) {
      GTEST_SKIP();
    }

    auto task = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTaskGetter)>(test_param)(
        this->GetTestInputData());

    EXPECT_FALSE(task->Validation());
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int size = std::get<0>(params);
    if (size > 0) {
      input_data_.resize(static_cast<size_t>(size));
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    static_cast<void>(output_data);
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

class KapanovaSIncreaseContrastEdgeTest : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int case_num = std::get<0>(params);

    switch (case_num) {
      case 1:
        input_data_ = {0};
        expected_output_ = {0};
        break;
      case 2:
        input_data_ = {255};
        expected_output_ = {255};
        break;
      case 3:
        input_data_ = {128, 128, 128};
        expected_output_ = {128, 128, 128};
        break;
      case 4:
        input_data_ = {0, 255, 0, 255};
        expected_output_ = {0, 255, 0, 255};
        break;
      case 5:
        input_data_ = {0, 64, 128, 192, 255};
        expected_output_ = {0, 64, 128, 192, 255};
        break;
      case 6:
        input_data_ = {100, 110, 120, 130};
        expected_output_ = {0, 85, 170, 255};
        break;
      case 7:
        input_data_ = {50, 100, 150, 200};
        expected_output_ = {0, 85, 170, 255};
        break;
      default:
        break;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_output_.size()) {
      return false;
    }

    for (size_t i = 0; i < output_data.size(); ++i) {
      if (std::abs(static_cast<int>(output_data[i]) - static_cast<int>(expected_output_[i])) > 1) {
        return false;
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

class KapanovaSIncreaseContrastRandomTest : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int size = std::get<0>(params);

    input_data_.resize(static_cast<size_t>(size));
    for (size_t i = 0; i < input_data_.size(); ++i) {
      input_data_[i] = static_cast<uint8_t>(50 + (i % 156));
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (input_data_.empty() || output_data.empty()) {
      return false;
    }
    if (input_data_.size() != output_data.size()) {
      return false;
    }

    auto min_in = *std::ranges::min_element(input_data_);
    auto max_in = *std::ranges::max_element(input_data_);
    auto min_out = *std::ranges::min_element(output_data);
    auto max_out = *std::ranges::max_element(output_data);

    if (min_in == max_in) {
      return std::ranges::all_of(output_data, [min_in](uint8_t v) { return v == min_in; });
    }

    if (min_out != 0 || max_out != 255) {
      return false;
    }

    for (size_t i = 0; i < input_data_.size(); ++i) {
      float scale = 255.0F / static_cast<float>(max_in - min_in);
      float normalized = static_cast<float>(input_data_[i] - min_in) * scale;
      auto expected = static_cast<uint8_t>(std::lround(normalized));

      if (std::abs(static_cast<int>(output_data[i]) - static_cast<int>(expected)) > 1) {
        return false;
      }
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(KapanovaSIncreaseContrastFuncTest, RunTest) {
  ExecuteTest(GetParam());
}

TEST_P(KapanovaSIncreaseContrastValidationTest, RunTest) {
  ExecuteTest(GetParam());
}

TEST_P(KapanovaSIncreaseContrastEdgeTest, RunTest) {
  ExecuteTest(GetParam());
}

TEST_P(KapanovaSIncreaseContrastRandomTest, RunTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParam = {std::make_tuple(4, "4x4"), std::make_tuple(8, "8x8"),
                                            std::make_tuple(16, "16x16"), std::make_tuple(32, "32x32"),
                                            std::make_tuple(64, "64x64")};

const std::array<TestType, 1> kTestValidationParam = {std::make_tuple(0, "empty")};

const std::array<TestType, 7> kTestEdgeParam = {std::make_tuple(1, "single_black"), std::make_tuple(2, "single_white"),
                                                std::make_tuple(3, "uniform_gray"), std::make_tuple(4, "bw_only"),
                                                std::make_tuple(5, "full_range"),   std::make_tuple(6, "narrow_range"),
                                                std::make_tuple(7, "linear_range")};

const std::array<TestType, 3> kTestRandomParam = {
    std::make_tuple(10, "random_small"), std::make_tuple(100, "random_medium"), std::make_tuple(1000, "random_large")};

const auto kTestTasksList = ppc::util::AddFuncTask<KapanovaSIncreaseContrastSEQ, InType>(
    kTestParam, PPC_SETTINGS_kapanova_s_increase_contrast_seq);

const auto kTestValidationTasksList = ppc::util::AddFuncTask<KapanovaSIncreaseContrastSEQ, InType>(
    kTestValidationParam, PPC_SETTINGS_kapanova_s_increase_contrast_seq);

const auto kTestEdgeTasksList = ppc::util::AddFuncTask<KapanovaSIncreaseContrastSEQ, InType>(
    kTestEdgeParam, PPC_SETTINGS_kapanova_s_increase_contrast_seq);

const auto kTestRandomTasksList = ppc::util::AddFuncTask<KapanovaSIncreaseContrastSEQ, InType>(
    kTestRandomParam, PPC_SETTINGS_kapanova_s_increase_contrast_seq);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kGtestValidationValues = ppc::util::ExpandToValues(kTestValidationTasksList);
const auto kGtestEdgeValues = ppc::util::ExpandToValues(kTestEdgeTasksList);
const auto kGtestRandomValues = ppc::util::ExpandToValues(kTestRandomTasksList);

const auto kTestName = KapanovaSIncreaseContrastFuncTest::PrintFuncTestName<KapanovaSIncreaseContrastFuncTest>;
const auto kValidationTestName =
    KapanovaSIncreaseContrastValidationTest::PrintFuncTestName<KapanovaSIncreaseContrastValidationTest>;
const auto kEdgeTestName = KapanovaSIncreaseContrastEdgeTest::PrintFuncTestName<KapanovaSIncreaseContrastEdgeTest>;
const auto kRandomTestName =
    KapanovaSIncreaseContrastRandomTest::PrintFuncTestName<KapanovaSIncreaseContrastRandomTest>;

INSTANTIATE_TEST_SUITE_P(IncreaseContrastTests, KapanovaSIncreaseContrastFuncTest, kGtestValues, kTestName);
INSTANTIATE_TEST_SUITE_P(ValidationTests, KapanovaSIncreaseContrastValidationTest, kGtestValidationValues,
                         kValidationTestName);
INSTANTIATE_TEST_SUITE_P(EdgeTests, KapanovaSIncreaseContrastEdgeTest, kGtestEdgeValues, kEdgeTestName);
INSTANTIATE_TEST_SUITE_P(RandomTests, KapanovaSIncreaseContrastRandomTest, kGtestRandomValues, kRandomTestName);

}  // namespace

}  // namespace kapanova_s_increase_contrast_seq
