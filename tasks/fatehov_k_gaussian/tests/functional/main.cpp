#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "fatehov_k_gaussian/common/include/common.hpp"
#include "fatehov_k_gaussian/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace fatehov_k_gaussian {

class FatehovKGaussianFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int size = std::get<0>(params);

    input_data_.image = Image(size, size, 3);
    for (size_t i = 0; i < input_data_.image.data.size(); i++) {
      input_data_.image.data[i] = static_cast<uint8_t>(i % 256);
    }
    input_data_.sigma = 1.0f;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.data.size() == input_data_.image.data.size();
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(FatehovKGaussianFuncTests, SequentialRun) {
  ExecuteTest(GetParam());
}

namespace {
const std::array<TestType, 3> kTestParam = {std::make_tuple(10, "10"), std::make_tuple(15, "15"),
                                            std::make_tuple(20, "20")};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<FatehovKGaussianSEQ, InType>(kTestParam, PPC_SETTINGS_fatehov_k_gaussian));

INSTANTIATE_TEST_SUITE_P(FatehovKTests, FatehovKGaussianFuncTests, ppc::util::ExpandToValues(kTestTasksList),
                         FatehovKGaussianFuncTests::PrintFuncTestName<FatehovKGaussianFuncTests>);
}  // namespace
}  // namespace fatehov_k_gaussian
