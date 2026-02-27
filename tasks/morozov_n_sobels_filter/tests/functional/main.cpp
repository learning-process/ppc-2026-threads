#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <tuple>

#include "morozov_n_sobels_filter/common/include/common.hpp"
#include "morozov_n_sobels_filter/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace morozov_n_sobels_filter {

class MorozovNRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return test_param;
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    std::string filename = params + ".txt";

    GetImageFromFile(filename);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if ((output_data.height != correct_image_.height) || (output_data.width != correct_image_.width) ||
        (output_data.pixels.size() != correct_image_.pixels.size())) {
      return false;
    }

    return output_data.pixels == correct_image_.pixels;
  }

  void GetImageFromFile(const std::string &path) {
    std::string abs_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_morozov_n_sobels_filter, path);
    std::ifstream file(abs_path);
    if (!file.is_open()) {
      throw std::runtime_error("failed to load image");
    }

    size_t height = 0;
    size_t width = 0;

    file >> height;
    file >> width;

    input_data_.height = height;
    input_data_.width = width;
    input_data_.pixels.resize(height * width);
    for (auto &pixel : input_data_.pixels) {
      int val = 0;
      file >> val;
      pixel = static_cast<uint8_t>(val);
    }

    std::string empty_line;
    std::getline(file, empty_line);

    correct_image_.height = height;
    correct_image_.width = width;
    correct_image_.pixels.resize(height * width);
    for (auto &pixel : correct_image_.pixels) {
      int val = 0;
      file >> val;
      pixel = static_cast<uint8_t>(val);
    }
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType correct_image_;
};

namespace {

TEST_P(MorozovNRunFuncTests, SobelOperator) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParam = {std::string("test_img_3x3"), std::string("test_img_5x5"),
                                            std::string("test_img_7x7"), std::string("test_img_9x9"),
                                            std::string("test_zero_9x9")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<MorozovNSobelsFilterSEQ, InType>(kTestParam, PPC_SETTINGS_morozov_n_sobels_filter));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = MorozovNRunFuncTests::PrintFuncTestName<MorozovNRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(SobelOperatorTests, MorozovNRunFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace morozov_n_sobels_filter
