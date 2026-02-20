#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "lukin_i_ench_contr_lin_hist/common/include/common.hpp"
#include "lukin_i_ench_contr_lin_hist/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace lukin_i_ench_contr_lin_hist {

class LukinIRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(test_param) + "x" + std::to_string(test_param) + "_size_image_test";
  }

 protected:
  void SetUp() override {
    TestType image_size = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int count = std::pow(image_size, 2);

    input_data_.resize(count);

    for (int i = 0; i < count; i++) {
      input_data_[i] = 80 + i % 81;  // [80,160] - как на обычных фото
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    unsigned char min = 255;
    unsigned char max = 0;

    for (const auto &elem : input_data_) {
      if (elem < min) {
        min = elem;
      }
      if (elem > max) {
        max = elem;
      }
    }

    float scale = 255.0f / (max - min);
    for (int i = 0; i < static_cast<int>(input_data_.size()); i++) {
      unsigned char expected_value = (input_data_[i] - min) * scale;
      if (output_data[i] != expected_value) {
        return false;
      }
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = {};
};

namespace {

TEST_P(LukinIRunFuncTestsThreads, LinearHist) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {128, 256, 512};  // размер изображения

const auto kTestTasksList =
    ppc::util::AddFuncTask<LukinITestTaskSEQ, InType>(kTestParam, PPC_SETTINGS_lukin_i_ench_contr_lin_hist);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = LukinIRunFuncTestsThreads::PrintFuncTestName<LukinIRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(PicLinearHist, LukinIRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace lukin_i_ench_contr_lin_hist
