#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include "kolotukhin_a_gaussian_blur/common/include/common.hpp"
#include "kolotukhin_a_gaussian_blur/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kolotukhin_a_gaussian_blur {

class KolotukhinAGaussinBlurePerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;
  std::string input_path_;

  void SetUp() override {
    const int rank = ppc::util::GetMPIRank();
    if (rank != 0) {
      input_data_ = {std::vector<std::uint8_t>{}, 0, 0};
      return;
    }
    input_path_ = ppc::util::GetAbsoluteTaskPath(PPC_ID_kolotukhin_a_gaussian_blur, "Frog_32x32.png");

    int width = -1;
    int height = -1;
    int channels = -1;
    unsigned char *raw = stbi_load(input_path_.c_str(), &width, &height, &channels, STBI_grey);
    if (raw == nullptr) {
      throw std::runtime_error("Load error: " + input_path_);
    }

    if (width <= 0 || height <= 0) {
      stbi_image_free(raw);
      throw std::runtime_error("Image has non-positive dimensions: " + input_path_);
    }

    const std::size_t img_size = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);

    std::vector<std::uint8_t> data(img_size);
    std::copy(raw, raw + img_size, data.begin());
    stbi_image_free(raw);

    input_data_ = {data, width, height};
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.size() == get<0>(input_data_).size();
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KolotukhinAGaussinBlurePerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KolotukhinAGaussinBlureSEQ>(PPC_SETTINGS_kolotukhin_a_gaussian_blur);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KolotukhinAGaussinBlurePerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KolotukhinAGaussinBlurePerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kolotukhin_a_gaussian_blur
