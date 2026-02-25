#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <new>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "kolotukhin_a_gaussian_blur/common/include/common.hpp"
#include "kolotukhin_a_gaussian_blur/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"
#include "util/include/util.hpp"

namespace kolotukhin_a_gaussian_blur {

class KolotukhinAGaussinBlurePerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    const int rank = ppc::util::GetMPIRank();
    if (rank != 0) {
      input_data_ = {std::vector<std::uint8_t>{}, 0, 0};
      return;
    }
    std::string input_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_kolotukhin_a_gaussian_blur, "test_image_1.jpg");

    int width = -1;
    int height = -1;
    int channels_in_file = -1;

    unsigned char *raw = stbi_load(input_path.c_str(), &width, &height, &channels_in_file, STBI_grey);

    if (raw == nullptr) {
      throw std::runtime_error("Load error: " + input_path);
    }

    if (width <= 0 || height <= 0) {
      stbi_image_free(raw);
      throw std::runtime_error("Image has non-positive dimension: " + input_path);
    }

    const std::size_t img_size = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);

    try {
      std::vector<std::uint8_t> data(img_size);
      std::copy(raw, raw + img_size, data.begin());
      input_data_ = {std::move(data), width, height};
    } catch (const std::bad_alloc &e) {
      stbi_image_free(raw);
      throw std::runtime_error("Failed to allocate memory for image: " + std::string(e.what()));
    }
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
