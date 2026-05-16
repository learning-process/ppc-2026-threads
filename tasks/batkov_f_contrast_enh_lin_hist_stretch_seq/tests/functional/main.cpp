#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <tuple>

#include "batkov_f_contrast_enh_lin_hist_stretch_seq/common/include/common.hpp"
#include "batkov_f_contrast_enh_lin_hist_stretch_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace batkov_f_contrast_enh_lin_hist_stretch_seq {

class BatkovFRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    const TestType &params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    auto size = static_cast<size_t>(std::get<0>(params));

    int base_intensity = 120;
    int contrast_range = 40;
    std::uniform_int_distribution<> dis(-contrast_range / 2, contrast_range / 2);

    input_data_.resize(size * size);
    for (size_t row = 0; row < size; ++row) {
      for (size_t col = 0; col < size; ++col) {
        int noise = dis(gen_);
        int value = base_intensity + noise;
        value = std::clamp(value, 80, 160);
        input_data_[(row * size) + col] = static_cast<uint8_t>(value);
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    auto [min_it, max_it] = std::ranges::minmax_element(output_data);
    return (*min_it == 0 && *max_it == 255) || (*min_it == *max_it);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  std::random_device rd_;
  std::mt19937 gen_{rd_()};
  InType input_data_;
};

namespace {

TEST_P(BatkovFRunFuncTestsThreads, ContrastEnhLinHistStretch) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 4> kTestParam = {std::make_tuple(100, "small_image"), std::make_tuple(500, "medium_image"),
                                            std::make_tuple(1000, "big_image"), std::make_tuple(2000, "large_image")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<BatkovFContrastEnhLinHistStretchSEQ, InType>(
    kTestParam, PPC_SETTINGS_batkov_f_contrast_enh_lin_hist_stretch_seq));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = BatkovFRunFuncTestsThreads::PrintFuncTestName<BatkovFRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(ContrastEnhLinHistStretch, BatkovFRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace batkov_f_contrast_enh_lin_hist_stretch_seq
