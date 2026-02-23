#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <tuple>

#include "gutyansky_a_img_contrast_incr/common/include/common.hpp"
#include "gutyansky_a_img_contrast_incr/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace gutyansky_a_img_contrast_incr {

class GutyanskyARunFuncTestsImgContrastIncr : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return test_param;
  }

 protected:
  void SetUp() override {
    std::string file_name =
        std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam()) + ".txt";
    std::string abs_path = ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_gutyansky_a_img_contrast_incr), file_name);

    std::ifstream ifs(abs_path);

    if (!ifs.is_open()) {
      throw std::runtime_error("Failed to open test file: " + file_name);
    }

    size_t num_of_elements = 0;
    ifs >> num_of_elements;

    input_data_.resize(num_of_elements);
    output_data_.resize(num_of_elements);

    for (auto &v : input_data_) {
      uint32_t value = 0;
      ifs >> value;
      v = static_cast<uint8_t>(value);
    }

    for (auto &v : output_data_) {
      uint32_t value = 0;
      ifs >> value;
      v = static_cast<uint8_t>(value);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (output_data_ == output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType output_data_;
};

namespace {

TEST_P(GutyanskyARunFuncTestsImgContrastIncr, ImgContrastIncr) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 10> kTestParam = {"test_1", "test_2", "test_3", "test_4", "test_5",
                                             "test_6", "test_7", "test_8", "test_9", "test_10"};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<GutyanskyAImgContrastIncrSEQ, InType>(
    kTestParam, PPC_SETTINGS_gutyansky_a_img_contrast_incr));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    GutyanskyARunFuncTestsImgContrastIncr::PrintFuncTestName<GutyanskyARunFuncTestsImgContrastIncr>;

INSTANTIATE_TEST_SUITE_P(ImgContrastIncrTests, GutyanskyARunFuncTestsImgContrastIncr, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace gutyansky_a_img_contrast_incr
