#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "otcheskov_s_contrast_lin_stretch/common/include/common.hpp"
#include "otcheskov_s_contrast_lin_stretch/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace otcheskov_s_contrast_lin_stretch {
namespace {

std::vector<uint8_t> CreateLowContrastImage(size_t size, uint8_t low = 100, uint8_t range = 50) {
  std::vector<uint8_t> image(size * size);
  for (size_t row = 0; row < size; ++row) {
    for (size_t col = 0; col < size; ++col) {
      uint8_t value = low + ((row + col) % range);
      image[(row * size) + col] = value;
    }
  }
  return image;
}

std::vector<uint8_t> LoadGrayscaleImage(const std::string &img_path) {
  int width = 0;
  int height = 0;
  int channels_in_file = 0;

  auto *data = stbi_load(img_path.c_str(), &width, &height, &channels_in_file, STBI_grey);
  if (data == nullptr) {
    throw std::runtime_error("Failed to load image '" + img_path + "': " + std::string(stbi_failure_reason()));
  }
  std::vector<uint8_t> img_data(data, data + static_cast<ptrdiff_t>(width * height));
  stbi_image_free(data);
  return img_data;
}

}  // namespace

class OtcheskovSContrastLinStretchValidationTestsThreads
    : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<0>(test_param) + "_" + std::to_string(std::get<1>(test_param));
  }

 protected:
  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.empty();
  }

  InType GetTestInputData() final {
    return input_data_;
  }

  void ExecuteTest(::ppc::util::FuncTestParam<InType, OutType, TestType> test_param) {
    const std::string &test_name =
        std::get<static_cast<std::size_t>(::ppc::util::GTestParamIndex::kNameTest)>(test_param);

    ValidateTestName(test_name);

    const auto test_env_scope = ppc::util::test::MakePerTestEnvForCurrentGTest(test_name);

    if (IsTestDisabled(test_name)) {
      GTEST_SKIP();
    }

    if (ShouldSkipNonMpiTask(test_name)) {
      std::cerr << "kALL and kMPI tasks are not under mpirun\n";
      GTEST_SKIP();
    }

    task_ =
        std::get<static_cast<std::size_t>(::ppc::util::GTestParamIndex::kTaskGetter)>(test_param)(GetTestInputData());
    const TestType &params = std::get<static_cast<std::size_t>(::ppc::util::GTestParamIndex::kTestParams)>(test_param);
    task_->GetInput() = std::vector<uint8_t>(std::get<1>(params));
    ExecuteTaskPipeline();
  }

  void ExecuteTaskPipeline() {
    EXPECT_FALSE(task_->Validation());
    task_->PreProcessing();
    task_->Run();
    task_->PostProcessing();
  }

 private:
  InType input_data_;
  ppc::task::TaskPtr<InType, OutType> task_;
};

class OtcheskovSContrastLinStretchFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<0>(test_param) + "_" + std::to_string(std::get<1>(test_param));
  }

 protected:
  void SetUp() override {
    const TestType &params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_img_.resize(std::get<1>(params));
    input_img_ = CreateLowContrastImage(std::get<1>(params));
  }

  bool CheckTestOutputData(OutType &output_data) final {
    auto [min_it, max_it] = std::ranges::minmax_element(output_data);
    return (*min_it == 0 && *max_it == 255) || (*min_it == *max_it);
  }

  InType GetTestInputData() final {
    return input_img_;
  }

 private:
  InType input_img_;
};

class OtcheskovSContrastLinStretchRealTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    std::string filename = std::get<0>(test_param);

    size_t dot_pos = filename.find_last_of('.');
    if (dot_pos != std::string::npos) {
      filename = filename.substr(0, dot_pos);
    }

    return filename;
  }

 protected:
  void SetUp() override {
    try {
      const TestType &params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
      const std::string &filename = std::get<0>(params);
      std::string abs_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_otcheskov_s_contrast_lin_stretch, filename);
      input_img_ = LoadGrayscaleImage(abs_path);
    } catch (const std::exception &e) {
      std::cerr << e.what() << '\n';
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    auto [min_it, max_it] = std::ranges::minmax_element(output_data);
    return (*min_it == 0 && *max_it == 255) || (*min_it == *max_it);
  }

  InType GetTestInputData() final {
    return input_img_;
  }

 private:
  InType input_img_;
};

namespace {

const std::array<TestType, 1> kTestValidParam = {{{"empty_data", 0}}};

const std::array<TestType, 5> kTestFuncParam = {
    {{"image_1x1", 1}, {"image_2x2", 2}, {"image_3x3", 3}, {"image_100x100", 100}, {"image_1000x1000", 1000}}};

const std::array<TestType, 1> kTestRealParam = {{{"grayimg.jpg", 0}}};

const auto kTestValidTasksList = std::tuple_cat(ppc::util::AddFuncTask<OtcheskovSContrastLinStretchSEQ, InType>(
    kTestValidParam, PPC_SETTINGS_otcheskov_s_contrast_lin_stretch));

const auto kTestFuncTasksList = std::tuple_cat(ppc::util::AddFuncTask<OtcheskovSContrastLinStretchSEQ, InType>(
    kTestFuncParam, PPC_SETTINGS_otcheskov_s_contrast_lin_stretch));

const auto kTestRealTasksList = std::tuple_cat(ppc::util::AddFuncTask<OtcheskovSContrastLinStretchSEQ, InType>(
    kTestRealParam, PPC_SETTINGS_otcheskov_s_contrast_lin_stretch));

const auto kGtestValidValues = ppc::util::ExpandToValues(kTestValidTasksList);
const auto kGtestFuncValues = ppc::util::ExpandToValues(kTestFuncTasksList);
const auto kGtestRealValues = ppc::util::ExpandToValues(kTestRealTasksList);

const auto kValidFuncTestName = OtcheskovSContrastLinStretchValidationTestsThreads::PrintFuncTestName<
    OtcheskovSContrastLinStretchValidationTestsThreads>;

const auto kFuncTestName =
    OtcheskovSContrastLinStretchFuncTestsThreads::PrintFuncTestName<OtcheskovSContrastLinStretchFuncTestsThreads>;

const auto kRealTestName =
    OtcheskovSContrastLinStretchRealTestsThreads::PrintFuncTestName<OtcheskovSContrastLinStretchRealTestsThreads>;

TEST_P(OtcheskovSContrastLinStretchValidationTestsThreads, ContrastLinStretchValidation) {
  ExecuteTest(GetParam());
}

TEST_P(OtcheskovSContrastLinStretchFuncTestsThreads, ContrastLinStretchFunc) {
  ExecuteTest(GetParam());
}

TEST_P(OtcheskovSContrastLinStretchRealTestsThreads, ContrastLinStretchReal) {
  ExecuteTest(GetParam());
}

INSTANTIATE_TEST_SUITE_P(ContrastLinStretchValidation, OtcheskovSContrastLinStretchValidationTestsThreads,
                         kGtestValidValues, kValidFuncTestName);

INSTANTIATE_TEST_SUITE_P(ContrastLinStretchFunc, OtcheskovSContrastLinStretchFuncTestsThreads, kGtestFuncValues,
                         kFuncTestName);

INSTANTIATE_TEST_SUITE_P(ContrastLinStretchReal, OtcheskovSContrastLinStretchRealTestsThreads, kGtestRealValues,
                         kRealTestName);

}  // namespace

}  // namespace otcheskov_s_contrast_lin_stretch
