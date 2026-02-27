#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>

#include "kazennova_a_sobel_operator/common/include/common.hpp"
#include "kazennova_a_sobel_operator/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kazennova_a_sobel_operator {

class KazennovaARunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    const TestType &params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    auto size = static_cast<size_t>(std::get<0>(params));

    // Создаём изображение с простым градиентом
    input_data_.resize(size * size);
    for (size_t i = 0; i < input_data_.size(); ++i) {
      input_data_[i] = static_cast<uint8_t>(i % 256);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // Проверка: выходное изображение не пустое и размер совпадает
    return !output_data.empty() && output_data.size() == input_data_.size();
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(KazennovaARunFuncTestsThreads, SobelTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(10, "small"), std::make_tuple(50, "medium"),
                                            std::make_tuple(100, "large")};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<SobelSeq, InType>(kTestParam, PPC_SETTINGS_kazennova_a_sobel_operator));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = KazennovaARunFuncTestsThreads::PrintFuncTestName<KazennovaARunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(SobelTests, KazennovaARunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kazennova_a_sobel_operator
